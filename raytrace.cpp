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
	TRANSMISSION,
	ALL
}Interaction;

//////////////////////////////////////////////////////////////////////
// Provides the framework for a raytracer.
////////////////////////////////////////////////////////////////////////
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
std::uniform_real_distribution<> myrandom(0 , 1.0f);
// Call myrandom(RNGen) to get a uniformly distributed random number in [0,1].
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

void Material::setTexture(const std::string path)
{
    //int width, height, n;
    //stbi_set_flip_vertically_on_load(true);
    //unsigned char* image = stbi_load(path.c_str(), &width, &height, &n, 0);

    //// Realtime code below:  This sends the texture in *image to the graphics card.
    //// The raytracer will not use this code (nor any features of OpenGL nor the graphics card).
    //glGenTextures(1, &texid);
    //glBindTexture(GL_TEXTURE_2D, texid);
    //glTexImage2D(GL_TEXTURE_2D, 0, n, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);

    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 100);
    //glGenerateMipmap(GL_TEXTURE_2D);
    //glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (int)GL_LINEAR);
    //glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (int)GL_LINEAR_MIPMAP_LINEAR);
    //glBindTexture(GL_TEXTURE_2D, 0);

    //stbi_image_free(image);
}

void Scene::Command(const std::vector<std::string>& strings, const std::vector<float>& f)
{
    if (strings.size() == 0) return;
    std::string c = strings[0];
    
    if (c == "screen") {
        // syntax: screen width height
      //  realtime->setScreen(int(f[1]),int(f[2]));
        width = int(f[1]);
        height = int(f[2]); }

    else if (c == "camera") {
        // syntax: camera x y z   ry   <orientation spec>
        // Eye position (x,y,z),  view orientation (qw qx qy qz),  frustum height ratio ry
      //  realtime->setCamera(Vector3f(f[1],f[2],f[3]), Orientation(5,strings,f), f[4]); 
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
        //realtime->setAmbient(Vector3f(f[1], f[2], f[3])); 
	}

    else if (c == "brdf")  {
        // syntax: brdf  r g b   r g b  alpha
        // later:  brdf  r g b   r g b  alpha  r g b ior
        // First rgb is Diffuse reflection, second is specular reflection.
        // third is beer's law transmission followed by index of refraction.
        // Creates a Material instance to be picked up by successive shapes
		//                           Kd                                Kr	              alpha
        currentMat = new Material(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7]); }

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
					 image[y * width + x] += tracer.TraceRay(r, lights, Tree);
				 }
			 }
		 }
		 printf(" %i",numberOfPasses);
		 numberOfPasses--;
	 } while (numberOfPasses > 0 && !Project1 );
	 printf("  at : %i ms\n",  (int)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count());
}


class Tracer
{
	void  SampleSphere(const Sphere* sphere, Intersection& result) {
		const float z = 2.0f * myrandom(RNGen) - 1.0f;
		const float r = sqrtf(1.0f - z * z);
		const float a = 2.0f * PI * myrandom(RNGen);
		result.n = Vector3f(r * cos(a), r * sin(a), z);
		result.p = sphere->center_ + result.n * sphere->radius;
	};

	////Phong BRDF
	Vector3f SampleBrdf(const Interaction & interaction, const Vector3f & wo, const Vector3f & normal, const float& alpha)
	{
		if (interaction == DIFFUSE) {
			return SampleLobe(normal, sqrtf(myrandom(RNGen)), 2 * PI * myrandom(RNGen));
		}
		else {// REFLECTION
			Vector3f m = SampleLobe(normal, pow(myrandom(RNGen), 1.0f / (alpha + 1.0f)), 2.0f * PI * myrandom(RNGen));
			return  2.0f * wo.dot(m) * m - wo;
		}
	}

	Vector3f  SampleLobe(Vector3f normal, float r1, float r2) {
		float const s = sqrtf(1 - r1 * r1);
		Vector3f const k = Vector3f(s * cosf(r2), s * sinf(r2), r1);
		Quaternionf const q = Quaternionf::FromTwoVectors(Vector3f::UnitZ(), normal);
		return q._transformVector(k);
	};

	float PdfBrdf( const  Vector3f & wo, Vector3f & normal, const  Vector3f & wi, Shape * shape) {
		float  Pd = fabsf(wi.dot(normal)) / PI;
		Vector3f m = (wo + wi).normalized();
		float Pr = D(normal, m, shape->material->alpha) * m.dot(normal) / (fabsf(wi.dot(m)));


		return Pd * shape->material->probabilityDiffuse +
			Pr * shape->material->probabilityReflection;
	};

	float PdfLight(const Intersection & result) {
		return 1.0f / result.s->get_area();
	};

	////Convert between angular measure and area measure
	float GeometryFactor(const Intersection & A, const Intersection & B) {
		const Vector3f  D = A.p - B.p;
		return fabsf(A.n.dot(D) * B.n.dot(D) / powf(D.dot(D), 2));
	};

	//Choose a uniformly distributed point on a sphere with center C and radius R
	Intersection  SampleLight(const std::vector<Shape*> & lights) {
		Intersection result;
		Sphere* sph = (Sphere*)lights[0];
		result.s = sph; //fills intersection.shape
		SampleSphere(sph, result);//fills intersection.Normal, intersection.Point
		return result;
	};

	Color F(const float& wiDOTm, const Color & Ks)//phong Fresnel
	{
		return Ks + (Color(1, 1, 1) - Ks) * (pow(1 - fabs(wiDOTm), 5));
	}
	float D(const Vector3f & normal, const Vector3f & m, const float& alpha)////Phong - distribution factor
	{
		return m.dot(normal) > 0.0f ? ((alpha + 2.0f) / (2.0f * PI)) * powf(m.dot(normal), alpha) : 0;
	}
	float G(const Vector3f & wo, const Vector3f & normal, const Vector3f & wi, const Vector3f & m, const float& alpha)//Phong geometric attenuation term
	{
		return G1(wi, m, normal, alpha)* G1(wo, m, normal, alpha);
	}
	float G1(const Vector3f & v, const  Vector3f & m, const  Vector3f & N, const  float& alpha) {
		float const vDOTn = v.dot(N);
		if (vDOTn > 1.0f)  return 1.0f;
		float const tanTheta = sqrtf(1.0f - (vDOTn * vDOTn)) / vDOTn;
		if (tanTheta == 0.0f) return 1.0f;

		float const a = sqrtf(((alpha / 2.0f) + 1.0f)) / fabsf(tanTheta);
		if (a < 1.6f) {
			float const  res = (3.535f * a + 2.181f * a * a) / (1.0f + 2.276f * a + 2.577f * a * a);
			return  v.dot(m) / vDOTn > 0.0f ? res : 0.0f;
		}
		else
			return v.dot(m) / vDOTn > 0.0f ? 1.0f : 0.0f;
	}

	Color EvalScattering(const Vector3f & wo, const Vector3f & normal, const Vector3f & wi, Shape * shape)
	{
		Color Ed = shape->material->Kd / PI;
		Vector3f m = (wo + wi).normalized();
		Color Er =
			D(normal, m, shape->material->alpha) *
			G(wo, normal, wi, m, shape->material->alpha) *
			F(wi.dot(m), shape->material->Ks) /
			(4.0f * fabsf(wi.dot(normal)) * fabsf(wo.dot(normal)));
		return fabsf(normal.dot(wi)) * (Ed + Er);
	};

public:
	Tracer() {};
	~Tracer() {};
	Intersection FindIntersection(const Ray & ray, const KdBVH<float, 3, Shape*> & tree) {
		Minimizer minimizer(ray);
		BVMinimize(tree, minimizer);
		return minimizer.minIntersection;
	};

	Color TraceRay(Ray & ray, std::vector<Shape*> & lights, const KdBVH<float, 3, Shape*> & tree) {
		Color C = Color(0.0f, 0.0f, 0.0f), W = Color(1.0f, 1.0f, 1.0f);
		Vector3f wi, wo = -ray.D;;//  
		float p;// PDF probability for color to happens
		Intersection P = FindIntersection(ray, tree);
		if (!P.s) return C; // no intersection
		if (P.s->isLight)  return P.s->material->Kd;// is light KD
		while (myrandom(RNGen) <= RUSSIAN_ROULETTE) {
			Vector3f normal = P.n;
			normal.normalize();
#ifdef EXPLICIT	// Explicit light connection
			Intersection L = SampleLight(lights);// Shadow Light
			p = PdfLight(L) / GeometryFactor(P, L);//Probability for that event to happen
			Vector3f wo = -ray.D;
			wi = L.p - P.p;
			P.p = P.p + wo * 0.0001f;
			Ray  I(P.p, wi);//Shadow Ray
			Intersection inters = FindIntersection(I, tree);
			if (p > 0 && inters.s && inters.p == P.p) {
				Color f = EvalScattering(wo, normal, wi, P.s);
				C += W * f / p * (Color)P.s->material->Kd;
			}
#endif
			Interaction myInteraction = myrandom(RNGen) < P.s->material->probabilityDiffuse ? DIFFUSE : REFLECTION;
			wi = SampleBrdf(myInteraction, wo, normal, P.s->material->alpha);// //Extend Path predict where the ray came from
			P.p += wi * 0.0001f; //fixing position outside shape. 
			Ray wiRay(P.p, wi);
			Intersection Q = FindIntersection(wiRay, tree);
			if (!Q.s) {
				C = Color(0.0f, 0.0f, 0.0f);
				break;
			}
			Color f = EvalScattering(wo, normal, wi, P.s);//
			p = PdfBrdf( wo, normal, wi, P.s) * RUSSIAN_ROULETTE; //how possible is that event
			if (p < epsilon) { //if there is not posibility
				C = Color(0.0f, 0.0f, 0.0f);
				break;
			}
			W *= f / p;
			if (Q.s->material->isLight()) {	//Implicit Light Connection
				C += W * (Color)Q.s->material->Kd;
				break;
			}
			P.Set(Q); //Prepare next possible ray rebound
			wo = -wi;
		}
		return C;
	};
};


//PDF probability density function
//PDFBRDF