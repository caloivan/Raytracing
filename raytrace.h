///////////////////////////////////////////////////////////////////////
// A framework for a raytracer.
////////////////////////////////////////////////////////////////////////
class Shape;
class Intersection;
const float PI = 3.14159265358979323846f;
const float Radians = PI / 180.0f;
const float epsilon = 0.000001f; //min distance offset for interseptions
////////////////////////////////////////////////////////////////////////
// Material: encapsulates a BRDF and communication with a shader.
////////////////////////////////////////////////////////////////////////
class Material
{
 public:
    Vector3f Kd, Ks;
    float alpha;
    unsigned int texid;
    virtual bool isLight() { return false; }

    Material()  : Kd(Vector3f(1.0, 0.5, 0.0)), Ks(Vector3f(1,1,1)), alpha(1.0), texid(0) {}
    Material(const Vector3f d, const Vector3f s, const float a) 
        : Kd(d), Ks(s), alpha(a), texid(0) {}
    Material(Material& o) { Kd=o.Kd;  Ks=o.Ks;  alpha=o.alpha;  texid=o.texid; }

    void setTexture(const std::string path);
    //virtual void apply(const unsigned int program);
};

////////////////////////////////////////////////////////////////////////
// Data structures for storing meshes -- mostly used for model files
// read in via ASSIMP.
//
// A MeshData holds two lists (stl::vector) one for vertices
// (VertexData: consisting of point, normal, texture, and tangent
// vectors), and one for triangles (TriData: consisting of three
// indices into the vertex array).
typedef Eigen::Matrix<unsigned int, 3, 1 > TriData;
    
class VertexData
{
 public:
    Vector3f pnt;
    Vector3f nrm;
    Vector2f tex;
    Vector3f tan;
    VertexData(const Vector3f& p, const Vector3f& n, const Vector2f& t, const Vector3f& a) 
        : pnt(p), nrm(n), tex(t), tan(a) 
    {}
};

struct MeshData
{
    std::vector<VertexData> vertices;
    std::vector<TriData> triangles;
    Material *mat;
};

class Light: public Material
{
public:

    Light(const Vector3f e) : Material() { Kd = e; }
    virtual bool isLight() { return true; }
    //virtual void apply(const unsigned int program);
};

class Realtime;
class Camera;

class Scene {
public:
    int width, height;
    //Realtime* realtime;         // Remove this (realtime stuff)
    Material* currentMat;
	Camera * camera_;
	std::vector<Shape*> shapes;
	std::vector<Shape*> lights;
    Scene();
	~Scene();
    void Finit();

    // The scene reader-parser will call the Command method with the
    // contents of each line in the scene file.
    void Command(const std::vector<std::string>& strings,
                 const std::vector<float>& f);

    // To read a model file into the scene via ASSIMP, call ReadAssimpFile.  
    void ReadAssimpFile(const std::string& path, const Matrix4f& M);

    // Once ReadAssimpFile parses the information from the model file,
    // it will call:
    void triangleMesh(MeshData* mesh);

    // The main program will call the TraceImage method to generate
    // and return the image.  This is the Ray Tracer!
    void TraceImage(Color* image, const int pass);

	Color ReturnColor(Intersection i, int option);
};

class Camera
{
public:
	Vector3f eye_;
	Quaternionf orientation_;
	float ry_;
	float rx_;
public:
	Camera(Vector3f eye, Quaternionf orientation, float ry, float width, float height) {
		eye_ = eye;
		orientation_ = orientation;
		ry_ = ry;
		rx_ = ry * (width / height);
	};
	~Camera();

};

class Ray {
public :
	Vector3f Q, D;
	Ray(Vector3f _Q, Vector3f _D) :Q(_Q), D(_D) {D = D.normalized();}
	Vector3f eval(float t) {return Q + t * D;}
};

class Intersection {
public:
	float t;//time of interception
	Shape* s;// shape interceted;
	Vector3f p;//point of interception
	Vector3f n;//normal of surface at interception
	float U, V;
public:
	Intersection() {
		t = INFINITY;//time of interception
		s = nullptr;// shape interceted;
	}

	void SetIntersection(float _time, Shape* _shape, Vector3f _P, Vector3f _N) {
		t = _time;
		p = _P;
		n = _N;
		s = _shape;
	}
	void SetIntersection(Intersection B) {
		t = B.t;
		s = B.s;
		p = B.p;
		n = B.n;
	}
};

class Slab {
public:
	float d0;		//offset 0
	float d1;		//offset 1
	Vector3f normal; //where  normal is pointing
	Slab(float _d0, float _d1, Vector3f _normal) :
		d0(_d0 > _d1 ? _d0 : _d1), d1(_d0 < _d1 ? _d0 : _d1), normal(_normal) {}
};

class Interval {
public:
	float  t0;// t min of interception
	float  t1;// t max of interception
public:
	Interval() : t0(0), t1(INFINITY) {}
	void Intersect(Ray _ray, const Slab& slab)// : forms interval by intersecting ray with slab, and intersects with this*
	{
		float nDotD = slab.normal.dot(_ray.D);
		float nDotQ = slab.normal.dot(_ray.Q);
		float s0 = slab.d0 + nDotQ;
		float s1 = slab.d1 + nDotQ;
		if (nDotD) { //ray direction and normal of plane have different directions. 
			t0 = -(s0) / (nDotD);
			t1 = -(s1) / (nDotD);
			if (t0 > t1) {
				float c = t0;
				t0 = t1;
				t1 = c;
			}
		}
		else {
			if (s0 * s1 < 0) {	// between the planes 
				t0 = 0;
				t1 = INFINITY;
			}
			else {     //outside the planes
				t0 = 1;
				t1 = 0;
			}
		}
	}
};

class Shape {
public:
	Shape():isLight(false) {  }
	virtual bool Intersect(Ray, Intersection&) = 0;
	virtual Bbox  returnBbox() = 0;
	virtual float get_area() = 0;
public:
	Material* material;
	Bbox boundingBox;
	bool isLight;
};

class Sphere : public Shape {
public:
	Vector3f center_;
	float radius;

public:
	Sphere() {}
	Sphere(Vector3f _center, float _radius, Material* _material) :center_(_center), radius(_radius)
	{
		material = _material;
		Vector3f radius3D(_radius, _radius, _radius);
		boundingBox = Bbox(_center - radius3D, _center + radius3D);
	}
	virtual bool Intersect(Ray ray, Intersection& _inter) {
		Vector3f q = ray.Q - center_;
		float qDOTd = q.dot(ray.D);
		float qDOTq = q.dot(q);
		float  discriminant = sqrtf(qDOTd * qDOTd - qDOTq + radius * radius);
		float tmin = -(qDOTd + discriminant),   tmax = -qDOTd + discriminant;
		float tempo = tmin > 0 ? tmin : tmax > 0 ? tmax :  0;
		if (tempo == 0) return false;
		Vector3f intersectionPoint = ray.eval(tempo);
		Vector3f normal = intersectionPoint -  center_;// (intersectionPoint – center_);
		_inter.SetIntersection(tempo, this, intersectionPoint, normal.normalized());// Object intercepted  
		return true;
	}
	Bbox  returnBbox() { return boundingBox; }
	float Sphere::get_area() { return 4 * PI* (radius * radius); }
};

class Box :public Shape{
public:
	Vector3f base;
	Vector3f diagonal;
public :
	Box() {}
	Box(Vector3f _base, Vector3f _diagonal, Material* _material) : base(_base), diagonal(_diagonal) {
		material = _material;
		boundingBox = Bbox(_base, _base + _diagonal);
	}
	virtual bool Intersect(Ray _ray, Intersection& _inter) {
		Slab slabX(-base[0], -base[0] - diagonal[0], Vector3f(1, 0, 0));//x
		Slab slabY(-base[1], -base[1] - diagonal[1], Vector3f(0, 1, 0));//y
		Slab slabZ(-base[2], -base[2] - diagonal[2], Vector3f(0, 0, 1));//z
		Interval intervalX, intervalY, intervalZ;
		intervalX.Intersect(_ray, slabX);//when detecting intersections, it returns intersection times. for both boundaries
		intervalY.Intersect(_ray, slabY);
		intervalZ.Intersect(_ray, slabZ);
		float tMi = fmax(intervalX.t0, fmax(intervalY.t0, intervalZ.t0));
		float tMa = fmin(intervalX.t1, fmin(intervalY.t1, intervalZ.t1));
		if (tMi > tMa) return false;// invalid behavior
		if (tMa <= epsilon) return false; //both behind the camera
		tMi = tMi > epsilon ? tMi : tMa; //minimum interception in front of the  camera
		Vector3f normal =
			tMi == intervalX.t0 ? Vector3f(1, 0, 0) :
			tMi == intervalX.t1 ? Vector3f(-1, 0, 0) :
			tMi == intervalY.t0 ? Vector3f(0, 1, 0) :
			tMi == intervalY.t1 ? Vector3f(0, -1, 0) :
			tMi == intervalZ.t0 ? Vector3f(0, 0, 1) :
			Vector3f(0, 0, -1);
		_inter.SetIntersection(tMi, this, _ray.eval(tMi), normal);
		return true;
	}
	Bbox  returnBbox() { return boundingBox; }
	virtual float get_area() { return 0.0f; };

};

class Cylinder: public Shape {
public :
	Vector3f base;// base;
	Vector3f axis;// axis;
	float radius;

public:
	Vector3f MinVec3(Vector3f& a, Vector3f& b)
	{
		return Vector3f(
			fmin(a[0], b[0]), fmin(a[1], b[1]), fmin(a[2], b[2]));
	}

	Vector3f MaxVec3(Vector3f& a, Vector3f& b)
	{
		return Vector3f(fmax(a[0], b[0]), fmax(a[1], b[1]), fmax(a[2], b[2]));
	}
	Cylinder(Vector3f _base, Vector3f _axis, float _radius, Material* _material) :
		base(_base), axis(_axis), radius(_radius) {
		material = _material;

		Vector3f v3r = Vector3f(radius, radius, radius);
		Vector3f lim1a = base - v3r;
		Vector3f lim1b = base + v3r;
		Vector3f lim2a = base + axis - v3r;
		Vector3f lim2b = base + axis + v3r;

		boundingBox = Bbox(
			MinVec3(lim2b, MinVec3(lim1b, MinVec3(lim1a, lim2a))),
			MaxVec3(lim1b, MaxVec3(lim2a, lim2b))
		);
	}

	virtual bool Intersect(Ray  _ray, Intersection& _inter)
	{
		Ray ray2 = _ray;
		float tLow[4];// array of  intersection times . to obtain the final slab
		float tHigh[4];
		// transform cylinder and ray to align the cylinder on Z-axis,  with endpoints at(0, 0, 0) and (0, 0, ‖A‖)
		Quaternionf q = Quaternionf::FromTwoVectors(axis, Vector3f::UnitZ());// cylinder's rotation
		//Change basis of ray to meet needs of cylinder , Transform ray Q + t D to 
		_ray.Q = q._transformVector(_ray.Q - base); //adjust ray to simulate cylinders rotation
		_ray.D = q._transformVector(_ray.D);
		float height = sqrtf( powf(axis[0] ,2) + powf(axis[1],2) + powf(axis[2] ,2));/*Cylinder's height from base to axis*/
		//intersect 3 intervals of the ray/*   
		
		/*  1.-    (0, infinity)*/
		tLow[1] = 0;
		tHigh[1] = INFINITY;

		/*   2.-   intersect with slabs on planes N(0,0,1) , d0 = 0 , d1= -||A|| */
		Slab slab(0, -height, Vector3f(0, 0, 1.0f));
		Interval interval;
		interval.Intersect(_ray, slab);
		tLow[2] = interval.t0;
		tHigh[2] = interval.t1;

		/*3.-   equation infinite cylinder  x^2+ y^2 = r^2, equation ray P(t) = Q + tD 
		mix equations (xE + t.xD)^2 + (yE + t.yD)^2 - r^2 = 0;
		solve for quadratic form  a.t^2 + b.t + c = 0,   t = ( -b +/- sqrt( b^2 - 4ac) ) ) /2a 
		( xE^2 + 2 * t.xD *xE + t^2 *xD^2 )   + *yE^2 + 2*yE*t.yD + t^2*t.yD^2 -r^2 = 0
		t^2 *(xD^2  + yD^2)        + t* ( 2.xD.xE + 2.yD.yE )    +   (xE^2 + yE^2 - r^2) = 0
					a                         b                              c  */
		float a = (_ray.D[0]) * (_ray.D[0]) + (_ray.D[1]) * (_ray.D[1]);
		float b = 2.0f * (_ray.D[0] * _ray.Q[0] + _ray.D[1] * _ray.Q[1]);
		float c = _ray.Q[0] * _ray.Q[0] + _ray.Q[1] * _ray.Q[1] - radius * radius;

		float discriminant = b * b - 4.0f * a * c;
		if (discriminant <= 0)   return false;// there weren't intersection points
		float sqrtDiscriminant = sqrtf(discriminant);
		tLow[3] = (-b - sqrtDiscriminant) / (2.0f * a);
		tHigh[3] = (-b + sqrtDiscriminant) / (2.0f * a);

		//step4 Calculate intersection of three intervals: [t0, t1] = [max(a0, b0, c0), min(a1, b1, c1)]
		tLow[0] = std::fmaxf(std::fmaxf(tLow[1], tLow[2]), tLow[3]);//define the maximum in three values
		tHigh[0] = std::fminf(std::fminf(tHigh[1], tHigh[2]), tHigh[3]);//define the minimum in three values
		
		if (tLow[0] >= tHigh[0]) return false; //step5 The “off the corner” case if limits arent cogruent
		if (tHigh[0] <= epsilon) return false; //object intercepts back the camera 

		float tIntersection = (tLow[0] > epsilon) ?tLow[0] : tHigh[0];//if t0 is positive and smaller than t1,if not, we know for sure t1 is positive and smaller than t0

		Vector3f intersectionPoint = _ray.eval(tIntersection); // calcule point of interception

		Vector3f normal = (tIntersection == tLow[2]) ? Vector3f::UnitZ() :// from up cylinder
			(tIntersection == tHigh[2]) ? -Vector3f::UnitZ() :// from down cylinder
			Vector3f(intersectionPoint[0], intersectionPoint[1], 0);//from side of the cylinder
		normal = normal.normalized();
		// normal = q.inverse.toRotationMatrix(normal); what next instruction does
		normal = q.conjugate()._transformVector(normal);// final normal of intersection. into world coordinates

		//Creating an Interception
		Vector3f InterceptionPoint = ray2.eval(tIntersection);
		_inter.SetIntersection(tIntersection, this, InterceptionPoint, normal);

		float theta = atan2(normal(1), normal(0));
		Vector2f uv(theta / (2.0f * PI), normal(2) / height);

		return  true;
	}


	Bbox  returnBbox() { return boundingBox; }
	virtual float get_area() { return 0.0f; };
};

class Triangle : public Shape {
public:
	Vector3f
		v0, v1, v2, n0, n1, n2;
public:
	Vector3f MinVec3(Vector3f& a, Vector3f& b) { return Vector3f(fmin(a[0], b[0]), fmin(a[1], b[1]), fmin(a[2], b[2])); }
	Vector3f MaxVec3(Vector3f& a, Vector3f& b) { return Vector3f(fmax(a[0], b[0]), fmax(a[1], b[1]), fmax(a[2], b[2])); }

	Triangle() {}
	Triangle(Vector3f _v0, Vector3f _v1, Vector3f _v2, Vector3f _n0, Vector3f _n1, Vector3f _n2, Material* _material)
		: v0(_v0), v1(_v1), v2(_v2), n0(_n0), n1(_n1), n2(_n2) {
		material = _material;
		boundingBox = Bbox(MinVec3(v0, MinVec3(v1, v2)), MaxVec3(v0, MaxVec3(v1, v2)));
	}

	virtual bool Intersect(Ray  _ray, Intersection& _inter) {
		Vector3f E1 = v1 - v0;
		Vector3f E2 = v2 - v0;
		Vector3f normal = E1.cross(E2).normalized();// SLIDES SAYS E2 X E1
		Vector3f p = _ray.D.cross(E2);
		float d = p.dot(E1);
		if (d == 0.0f) return false;// paralel to triangle
		float dInv = 1.0f / d;
		Vector3f S = _ray.Q - v0;
		float u = p.dot(S) * dInv;
		if (u < 0.0f || u > 1.0f) return false;//  intersects triangle, but outside E2 edge

		Vector3f q = S.cross(E1);
		float v = _ray.D.dot(q) * dInv;
		if (v < 0.0f || (u + v) > 1.0f) return false; //ray intersects plane but outside other edges

		float t = E2.dot(q) * dInv;
		if (t < epsilon) return false; //  intersect behind eeye

		_inter.SetIntersection(t, this, _ray.eval(t), normal);
		return  true;
	}
	Bbox  returnBbox() { return boundingBox; }
	virtual float get_area() { return 0.0f; };
};