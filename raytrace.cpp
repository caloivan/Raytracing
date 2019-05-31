int optionDraw = 4; //1 deep 2 normal 3 material  4 KD
int kdtree = 1;
int numberOfPasses = 50;  // if passes = 1 then raytracing, else path tracing
bool Boxes = true;
bool Cylinders = true;
bool Spheres = true;

bool Project1 = false;
#define EXPLICIT

typedef enum {
	DIFFUSE,
	REFLECTION,
	TRANSMISSION
}Interaction;

#include <vector>
#include <windows.h>
#include <cstdlib>
#include <limits>
#include <crtdbg.h>
#include "geom.h"
#include "raytrace.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
#include <chrono>// measure time
std::mt19937_64 RNGen,RNGen2;// A good quality *thread-safe* Mersenne Twister random number generator.
std::uniform_real_distribution<> myrandom(0 , 1.0f);// Call myrandom(RNGen) to get a uniformly distributed random number in [0,1].
///************************** KDTREE SECTION ***************************************///
//These includes are  written on geom.h
//#include <Eigen/StdVector> // For vectors, matrices (2d,3d,4d) and quaternions in f and d precision.
//#include <Eigen_unsupported/Eigen/BVH> // For KdBVH
//typedef AlignedBox<float,3> Bbox;
class BoxVolume : public Shape {
public:
	Vector3f base;// base;
	Vector3f diagonal;// diagonal;
public:
	BoxVolume(Vector3f _b0, Vector3f _b1) :base(_b0), diagonal(_b1) {}
	virtual bool Intersect(Ray _ray, Intersection & _inter) {
		Slab slabX(-base[0], -base[0] - diagonal[0], Vector3f(1, 0, 0));
		Slab slabY(-base[1], -base[1] - diagonal[1], Vector3f(0, 1, 0));
		Slab slabZ(-base[2], -base[2] - diagonal[2], Vector3f(0, 0, 1));

		Interval intervalX, intervalY, intervalZ;
		intervalX.Intersect(_ray, slabX);
		intervalY.Intersect(_ray, slabY);
		intervalZ.Intersect(_ray, slabZ);

		float tMi = fmax(intervalX.t0, fmax(intervalY.t0, intervalZ.t0));
		float tMa = fmin(intervalX.t1, fmin(intervalY.t1, intervalZ.t1));

		if (tMi > tMa || tMa <= epsilon) return false;// invalid behavior  //both behind ray

		tMi = tMi > epsilon ? tMi : 0; //minimum interception in front of camera, if ray starts inside box, return 0
		_inter.Set(tMi, this, _ray.eval(tMi), Vector3f(0, 0, 0));// normal doesnt matters, cause here is a rough intersection

		return true;
	}
	Bbox  returnBbox() { return boundingBox; }
	virtual float get_area() { return 0.0f; };
};

class Minimizer
{
public:
	typedef float Scalar; // KdBVH needs Minimizer::Scalar defined
	Ray ray;
	Intersection minIntersection;
public:
	Minimizer(const Ray& r) : ray(r) { }
	float minimumOnObject(Shape* obj) {
		Intersection  intersection;
		if (obj->Intersect(ray, intersection)) {
			minIntersection = intersection.t < minIntersection.t ? intersection : minIntersection;
			return intersection.t;
		}
		return INFINITY;
	}
	float minimumOnVolume(const Bbox& box) {
		BoxVolume volume((box.min)(), (box.max)() - (box.min)());
		Intersection i;
		return volume.Intersect(ray, i) ? i.t : INFINITY;//if volume intersects ray, return time , else return infinity
	}
};

Bbox bounding_box(Shape* shape) { return shape->returnBbox(); }
///************************** END SECTION ***************************************///

class Tracer
{
	Intersection FindIntersection(const Ray& ray, const KdBVH<float, 3, Shape*>& tree ) {
		Minimizer minimizer(ray);
		BVMinimize(tree, minimizer);
		return minimizer.minIntersection;
	};

public: 
Color Tracer::TraceRay(Ray& ray, std::vector<Shape*>& lights, KdBVH<float, 3, Shape*>& tree) {
	Color C = Color(0.0f, 0.0f, 0.0f);//cOLOR
	Color W = Color(1.0f, 1.0f, 1.0f);// WEIGHT PONDERATOR
	float p, q, Wmis;
	bool objectMoved = false; //reference to move 1 object
	Vector3f offset = Vector3f(0, 0, 0);
	Vector3f initialPosition = Vector3f(0, 0, 0);
	//INITIAL RAY
	Intersection P = FindIntersection(ray, tree);
	if (!P.s) return C; // no intersection
	if (P.s->isLight)  return P.s->material->Kd;// is light KD
	Vector3f wi, wo = -ray.D;
	while (myrandom(RNGen) <= RUSSIAN_ROULETTE) {
		Vector3f normal = P.n.normalized();
#ifdef EXPLICIT
		Intersection L = SampleLight(lights);// Randomly choose a light 
		p = PdfLight(L) / GeometryFactor(P, L);//  Probability to hit light ( angular measure )
		q = PdfBrdf( wo, normal, wi, P.s) * RUSSIAN_ROULETTE;//probability of diffuse + reflection + refraction
		Wmis =  p * p / (p * p + q * q);
		wi = (L.p - P.p).normalized();
		Ray shadowRay((P.p - wo * 0.01f), wi);//Ray goes explicit to light
		Intersection inters = FindIntersection(shadowRay, tree);
		if (p > 0 && inters.s && inters.p == L.p) {
			Color f = EvalScattering(  wo, normal, wi, P.s, P.t);
			C += W * (f / p) * Wmis * (Color)L.s->material->Kd;
		}
#endif
		float rand = myrandom(RNGen);
		Interaction choice =
			rand < P.s->material->probabilityDiffuse ? DIFFUSE :
			rand < P.s->material->probabilityDiffuse + P.s->material->probabilityReflection ? REFLECTION : TRANSMISSION;
		wi = SampleBrdf(choice, normal, wo, P.s);// calculates rebound, based on diffuse, reflection, refraction
		Ray wiRay((P.p + wi * 0.0001f), wi);
		Intersection Q = FindIntersection(wiRay, tree);
		if (!Q.s)   break;

		Color f = EvalScattering( wo, normal, wi, P.s, P.t);
		p = PdfBrdf( wo, normal, wi, P.s) * RUSSIAN_ROULETTE;
		if (p < epsilon)   break;
		W *= f / p;

		//IMPLICIT LIGHT CONNECTION
		if (Q.s->material->isLight()) {
			q = PdfLight(Q) / GeometryFactor(P, Q);//Probability the implicit light could be chosen explicitly
			Wmis = 1;// p * p / (p * p + q * q);
			C += W * Wmis * (Color)Q.s->material->Kd;
			break;
		}
		P.Set(Q);
		wo = -wi;
	}
	return  C;
}

Color Tracer::EvalScattering( Vector3f & wo, Vector3f& normal, Vector3f & wi, Shape * shape, float& t) {
	    Vector3f Kd = shape->material->Kd;
	    Vector3f Ks = shape->material->Ks;
		Vector3f Kt = shape->material->Kt;
		float alpha = shape->material->alpha;
		float ior = shape->material->ior;
		Color ColorDiffuse = (Color)Kd / PI;
		
		float woDOTn = fabsf(wo.dot(normal));
		float wiDOTn = fabsf(wi.dot(normal));
		Vector3f m = (wo + wi).normalized();
		float wiDOTm = wi.dot(m);
		Color ColorReflection = F(wiDOTm, Ks) * D(normal, m, alpha) * G(wo, normal, wi, m,  alpha) / 
											 (4.0f * woDOTn * wiDOTn);

		float ni = wo.dot(normal) < 0.0f?ior: 1.0f;
		float no = wo.dot(normal) < 0.0f?1.0f:ior;
		float n = ni / no;
		Vector3f attuenation_ = wo.dot(normal) < 0.0f ?// BEERS LAW FOR ATENUATION
			Vector3f(powf(2.71f, t * log(Kt.x())), powf(2.71f, t * log(Kt.y())), powf(2.71f, t * log(Kt.z()))):
			Vector3f(1, 1, 1);
		 m = -(wo * ni + wi * no).normalized();
		float woDOTm = wo.dot(m);
		 wiDOTm = wi.dot(m);
		Color ColorTransmission;
		float radicand = 1 - (n * n) * (1 - (woDOTm * woDOTm));
		if (radicand < 0) {//INTERNAL REFLECTION
			ColorTransmission = F(wiDOTm, Ks) * D(normal, m, alpha) * G(wo, normal, wi, m, alpha)  / 
											 	(4.0f * woDOTn * wiDOTn);
		}
		else {
			ColorTransmission = fabsf(wiDOTm) * fabsf(woDOTm) * no * no      *      (Color(1, 1, 1) - F(wiDOTm, Ks)) * D(normal, m, alpha) * G(wo, normal, wi, m, alpha) /
				                  (powf(no * wiDOTm + ni * woDOTm, 2)        *	        woDOTn * wiDOTn);
		}
		ColorTransmission *= (Color)attuenation_;
	return fabsf(normal.dot(wi)) * (ColorDiffuse + ColorReflection + ColorTransmission);
}

//Returns the Wi. acording to normal. wo and  surface specig
Vector3f  Tracer::SampleBrdf(const Interaction& choice, Vector3f& normal, Vector3f& wo, Shape* shape) {
	Vector3f wi;
	float alpha = shape->material->alpha;
	if (choice == DIFFUSE) {
		wi = SampleLobe(normal, sqrtf(myrandom(RNGen)), 2 * PI * (myrandom(RNGen)));
	}
	else if (choice == REFLECTION) {
		Vector3f m = SampleLobe(normal, pow(myrandom(RNGen), 1.0f / (alpha + 1.0f)), 2 * PI * (myrandom(RNGen)));
		wi = 2.0f * wo.dot(m) * m - wo;
	}
	else if (choice == TRANSMISSION) {
		Vector3f m = SampleLobe(normal, pow(myrandom(RNGen), 1.0f / (alpha + 1.0f)), 2 * PI * (myrandom(RNGen)));
		float woDOTn = wo.dot(normal);
		float n =  woDOTn < 0 ? shape->material->ior :1/ shape->material->ior;
		float woDOTm = wo.dot(m);
		float radicand = 1.0f - (n * n) * (1 - (woDOTm * woDOTm));
		if (radicand < 0)
			wi = 2.0f * woDOTm * m - wo;
		else
			wi = (n * woDOTm - (woDOTn >= 0 ? 1 : -1) * sqrtf(radicand)) * m - n * wo;
	}
	return wi.normalized();
}

Vector3f Tracer::SampleLobe(Vector3f normal, float r1, float r2) {
	float s = sqrtf(1 - (r1 * r1));
	Quaternionf q = Quaternionf::FromTwoVectors(Vector3f::UnitZ(), normal);
	return q._transformVector(Vector3f(s * cosf(r2), s * sinf(r2), r1));
}

float Tracer::PdfBrdf(  Vector3f & wo, Vector3f & normal, Vector3f & wi, Shape * shape) {
	float alpha = shape->material->alpha;
	float pd = shape->material->probabilityDiffuse;
	float pr = shape->material->probabilityReflection;
	float pt = shape->material->probabilityTransmission;
	float Pd, Pt, Pr;
	Pd = fabsf(wi.dot(normal)) / PI;
	Vector3f m = (wo + wi).normalized();
	Pr =  D(normal, m, alpha) * (fabsf(m.dot(normal))) / (4.0f * fabsf(wi.dot(m)));
	
	float ior = shape->material->ior;
	float woDOTn = wo.dot(normal);
	float ni = woDOTn < 0 ? ior : 1.0f;
	float no = woDOTn < 0 ? 1 : ior;
	float n =  ni / no;
	m = -(wo * ni + wi * no).normalized();
	float woDOTm = wo.dot(m);
	float radicand = 1 - (n * n) * (1 - (woDOTm * woDOTm));

	m = radicand < 0.0f?(wo + wi).normalized() : m;
	float deno = radicand < 0.0f ? 1:(no * (wi.dot(m))) + (ni * (woDOTm));

	Pt = radicand < 0.0f? //total internal reflection
			D(normal, m, alpha) * fabsf(m.dot(normal)) / 
				(4.0f * fabsf(wi.dot(m))):

		    D(normal, m, alpha) * fabsf(m.dot(normal)) * no * no * fabsf(wi.dot(m))  /
				( deno * deno);
	return pd * Pd + pr * Pr + pt * Pt;
}

float Tracer::GeometryFactor(Intersection & A, Intersection & B) {
	const Vector3f  D = A.p - B.p;
			return fabsf(A.n.dot(D) * B.n.dot(D) / 
				        powf(D.dot(D), 2));
}

Intersection Tracer::SampleLight(std::vector<Shape*> & lights) {
	int index = myrandom(RNGen) * lights.size();
	return SampleSphere((Sphere*)lights[index]);
}

Intersection Tracer::SampleSphere(Sphere * sph) {
	Intersection result;
	float z = 2 * myrandom(RNGen) - 1.0f;
	float r = sqrtf(1 - z * z);
	float a = 2 * PI * myrandom(RNGen);
	result.n = Vector3f(r * cos(a), r * sin(a), z);
	result.p = sph->center_ + result.n * sph->radius;
	result.s = sph;
	return result;
}

float Tracer::PdfLight(Intersection & result) {
	return       1.0f / 
		    result.s->get_area();
}

Color Tracer::F(const float& wiDOTm, Vector3f Ks){// Shape* shp) {//TODO since  when wiDOTm equal to LDOTH
	return Ks + (Vector3f(1, 1, 1) - Ks) * pow(1 - fabs(wiDOTm), 5);
}

float Tracer::D(Vector3f & normal, Vector3f & m, float& alpha) {
	return ((m.dot(normal)) > 0 ? 1 : 0) * ((alpha + 2.0f) / (2.0f * PI)) * pow((m.dot(normal)), alpha);
}

float Tracer::G(Vector3f & wo, Vector3f& normal, Vector3f & wi, Vector3f & m,  float& alpha) {
	return G1(wi, m, normal, alpha)* G1(wo, m, normal, alpha);
}

float Tracer::G1(Vector3f & v, Vector3f & m, Vector3f & normal, float& alpha) {
	float vDOTn = v.dot(normal);
	if (vDOTn > 1.0f)   return 1.0f;

	float tanTheta = sqrtf(1.0f - (vDOTn * vDOTn)) / vDOTn;
	if (tanTheta == 0.0f)   return 1.0f;

	float a = sqrtf((alpha / 2.0f) + 1.0f) / fabs(tanTheta);

	if (a < 1.6f)
		return ((v.dot(m) / vDOTn) > 0 ? 1.0f : 0) * (3.535f * a + 2.181f * a * a) / (1.0f + 2.276f * a + 2.577f * a * a);
	else
		return ((v.dot(m) / vDOTn) > 0 ? 1.0f : 0);
}
};//Tracer Class

Scene::Scene() { }
void Scene::Finit(){}
void Scene::triangleMesh(MeshData* mesh)
{
	for (TriData& tri : mesh->triangles) {
		Triangle* pTri = new Triangle(
			mesh->vertices[tri[0]].pnt,
			mesh->vertices[tri[1]].pnt,
			mesh->vertices[tri[2]].pnt,
			Vector3f(0.0f, 0.0f, 0.0f),// TODO: need  for later
			Vector3f(0.0f, 0.0f, 0.0f),
			Vector3f(0.0f, 0.0f, 0.0f),
			currentMat
		);
		shapes.push_back(pTri);
	}
}
Quaternionf Orientation(unsigned int i,  const std::vector<std::string>& strings, const std::vector<float>& f)
{
    Quaternionf q(1,0,0,0); // Unit quaternion
    while (i<strings.size()) {
        std::string c = strings[i++];
        if (c == "x")  
            q *= angleAxis(f[i++]*Radians, Vector3f::UnitX());
        else if (c == "y")  
            q *= angleAxis(f[i++]*Radians, Vector3f::UnitY());
        else if (c == "z")  
            q *= angleAxis(f[i++]*Radians, Vector3f::UnitZ());
        else if (c == "q")  {
            q *= Quaternionf(f[i+0], f[i+1], f[i+2], f[i+3]);
            i+=4; }
        else if (c == "a")  {
            q *= angleAxis(f[i+0]*Radians, Vector3f(f[i+1], f[i+2], f[i+3]).normalized());
            i+=4; } }
    return q;
}

void Scene::Command(const std::vector<std::string>& strings, const std::vector<float>& f)
{
    if (strings.size() == 0) return;
    std::string c = strings[0];
    
    if (c == "screen") {
        // syntax: screen width height
        width = int(f[1]);
        height = int(f[2]); }

    else if (c == "camera") {
        // syntax: camera x y z   ry   <orientation spec>
        // Eye position (x,y,z),  view orientation (qw qx qy qz),  frustum height ratio ry
		camera_ = new Camera(
			Vector3f(f[1], f[2], f[3]),//Eye position (x,y,z)
			Orientation(5, strings, f),//view orientation (qw qx qy qz),
			f[4], // frustum height = 0.2f
			(float)width,
			(float)height
		);
	}

    else if (c == "ambient") {
        // syntax: ambient r g b
        // Sets the ambient color.  Note: This parameter is temporary.
        // It will be ignored once your raytracer becomes capable of
        // accurately *calculating* the true ambient light.
	}

	else if (c == "brdf") {
		// syntax: brdf  r g b   r g b  alpha
		// later:  brdf  r g b   r g b  alpha  r g b ior
		// First rgb is Diffuse reflection, second is specular reflection.
		// third is beer's law transmission followed by index of refraction.
		// Creates a Material instance to be picked up by successive shapes
		if (f.size() == 12)
			//                                       Kd                             Kr	                          Ks              alpha
			currentMat = new Material(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7], Vector3f(f[8], f[9], f[10]), f[11]);
		else
			//                           Kd                                Kr	              alpha
			currentMat = new Material(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7]);
	}
    else if (c == "light") {
        // syntax: light  r g b   
        // The rgb is the emission of the light
        // Creates a Material instance to be picked up by successive shapes
        currentMat = new Light(Vector3f(f[1], f[2], f[3])); }
   
    else if (c == "sphere" && Spheres) {
        // syntax: sphere x y z   r
        // Creates a Shape instance for a sphere defined by a center and radius
		Shape* mySPhere = new Sphere(Vector3f(f[1], f[2], f[3]),    f[4],    currentMat);
		shapes.push_back(mySPhere);
		if (currentMat->isLight()) {
			mySPhere->isLight = true;
			lights.push_back(mySPhere);
		}
	}

    else if (c == "box" && Boxes) {
        // syntax: box bx by bz   dx dy dz
        // Creates a Shape instance for a box defined by a corner point and diagonal vector
		Shape* myBox = new Box(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), currentMat);
		shapes.push_back(myBox);
	}

    else if (c == "cylinder" && Cylinders) {
        // syntax: cylinder bx by bz   ax ay az  r
        // Creates a Shape instance for a cylinder defined by a base point, axis vector, and radius
		Shape* myCylinder = new Cylinder(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7], currentMat);
		shapes.push_back(myCylinder);
	}

    else if (c == "mesh") {
        // syntax: mesh   filename   tx ty tz   s   <orientation>
        // Creates many Shape instances (one per triangle) by reading
        // model(s) from filename. All triangles are rotated by a
        // quaternion (qw qx qy qz), uniformly scaled by s, and
        // translated by (tx ty tz) .
        Matrix4f modelTr = translate(Vector3f(f[2],f[3],f[4]))
                          *scale(Vector3f(f[5],f[5],f[5]))
                          *toMat4(Orientation(6,strings,f));
        ReadAssimpFile(strings[1], modelTr);  }

    else {
        fprintf(stderr, "\n*********************************************\n");
        fprintf(stderr, "* Unknown command: %s\n", c.c_str());
        fprintf(stderr, "*********************************************\n\n");
    }
}

Color Scene::ReturnColor(Intersection i, int option) {
	switch (option) {
	case 1: {
		float myDeepth = (i.t - 5) * 0.25f;
		return Color(myDeepth, myDeepth, myDeepth);// DEEPTH 
	}break;
	case 2: {	
		return i.n.cwiseAbs();// NORMALS 
	}break;
	case 3: {
		return i.s->material->Kd;// MATERIAL
	}break;
	case 4: {
		Vector3f L = static_cast<Sphere*>(lights.front())->center_ - i.p;
		return i.n.dot(L)* i.s->material->Kd / PI;// DIFFUSE
	}break;
	default: break;
	}
	return Color(1.0f, 1.0f, 0.0f);
}

void Scene::TraceImage(Color* image, const int pass)
{
	 Vector3f const X =         camera_->rx_ * camera_->orientation_._transformVector(Vector3f::UnitX());
	 Vector3f const Y =         camera_->ry_ * camera_->orientation_._transformVector(Vector3f::UnitY());
	 Vector3f const Z = -1.0f *                camera_->orientation_._transformVector(Vector3f::UnitZ());

	 KdBVH<float, 3, Shape*> Tree(shapes.begin(), shapes.end());
	 auto const start = std::chrono::steady_clock::now();
	 printf("Pass \n");
	 do {
		#pragma omp parallel for schedule(dynamic, 1) // Magic: Multi-thread y loop
		 for (int y = 0; y < height; y++) {
			 for (int x = 0; x < width; x++) {
				 float const dx = 2.0f * ((x + (float)(myrandom(RNGen))) / (float)width) - 1.0f;
				 float const dy = 2.0f * ((y + (float)(myrandom(RNGen))) / (float)height) - 1.0f;
				 Ray r(camera_->eye_, dx * X + dy * Y + Z);
				 if (Project1) {
					 if (kdtree) {
						 Minimizer minimizer(r);
						 BVMinimize(Tree, minimizer);
						 if (minimizer.minIntersection.t != INFINITY)
							 image[y * width + x] = ReturnColor(minimizer.minIntersection, optionDraw);
					 }
					 else {
						 float minTime = INFINITY;
						 for (unsigned int i = 0; i < shapes.size(); ++i) {
							 Intersection intersection;
							 if (shapes[i]->Intersect(r, intersection)) {
								 if (intersection.t < minTime) {// compare 
									 minTime = intersection.t;
									 image[y * width + x] = ReturnColor(intersection, optionDraw);
								 }
							 }
						 }
					 }
				 }
				 else {
					 Tracer tracer;
					 image[y * width + x] += tracer.TraceRay(r, lights,  Tree);
				 }
			 }
		 }
		 printf(" %i",numberOfPasses);
		 numberOfPasses--;
	 } while (numberOfPasses > 0 && !Project1 );
	 printf("  at : %i ms\n",  (int)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count());
}

//PDF probability density function
//PDFBRDF