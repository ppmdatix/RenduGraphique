#include <iostream>

using namespace std;

#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>
#include <random>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "algorithm"

#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>


//Let's define what is a vector.
class Vector{
public:
    Vector(double x = 0, double y = 0, double z = 0): x(x), y(y), z(z) {};
    double norm2() {return x*x +y*y + z*z;}
    void normalize() {
        double n = sqrt(norm2());
        x/= n;
        y/=n;
        z/=n;
    }
    double x,y,z;

};
//Let's define oerations between vectors and or scalars
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a.x + b.x, a.y +b.y, a.z +b.z);
}

Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a.x - b.x, a.y -b.y, a.z -b.z);
 
}

Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a.x * b.x, a.y *b.y, a.z *b.z);
 
}

Vector operator-(const Vector& a) {
    return Vector(-a.x, -a.y, -a.z);
}

Vector operator*(const double b,const Vector& a) {
    return Vector(a.x * b, a.y *b, a.z *b);
}

Vector operator/(const Vector& a, const double b) {
    return Vector(a.x / b, a.y /b, a.z /b);
}

double dot(const Vector& a, const Vector& b) {
    return a.x * b.x+ a.y * b.y + a.z *b.z;
};

Vector cross(const Vector& a, const Vector& b) {
	return Vector(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
}

// A ray is an origin and a direction
class Ray {
public:
	Ray(const Vector& C , const Vector& u): C(C), u(u) {};
	
	Vector C;
	Vector u;
}	;


class Object{
public:
	Object() {};
	
	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t) = 0;
	
	Vector albedo;
    bool mirror;
    bool transparent;
    double nsphere;
    double emissivite;
};

class Sphere : public Object{
public:

    Sphere(const Vector& O, double R, Vector albedo, bool mirror, bool transparent, double nsphere, double emissivite): O(O), R(R){
		this->albedo = albedo;
		this->mirror = mirror;
		this->transparent = transparent;
		this->nsphere = nsphere;
		this->emissivite = emissivite;};
    
    virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t){
		double a = 1;
		double b = 2 * dot(r.u, r.C-O);
		double c = (r.C - O).norm2() - R*R;
		double delta  = b*b - 4 * a* c;
		P = Vector(0,0,0);
		if (delta>=0){
			double sqrtdelta = sqrt(delta);
			double t1 = (-b - sqrtdelta)/(2.*a);
			double t2 = (-b + sqrtdelta)/(2.*a);
			if (t1 > 0) {
				t=t1;
				P = t1*r.u + r.C;
			} 
			else {
				if (t2 > 0) {
					t=t2;
					P = t2*r.u + r.C;
				} 
				else {return false;}
			}
		}
		N = P- O;
		N.normalize();
		return (delta>=0);
	}
	
    Vector O;
    double R;

};
class Triangle : public Object{
public:

    Triangle(Vector A, Vector B, Vector C, Vector albedo, bool mirror, bool transparent, double nsphere, double emissivite): A(A), B(B), C(C){
		this->albedo = albedo;
		this->mirror = mirror;
		this->transparent = transparent;
		this->nsphere = nsphere;
		this->emissivite = emissivite;
		}
    
    virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t){
		//cout << "inter" << endl;
		N = -1.0 *  cross(B-A,B-C);
		N.normalize();
		if (dot(N,r.u) > 0){N = - N;}
		if (dot(r.u, N) == 0){return false;}
		t = dot(A-r.C, N) / dot(r.u, N);
		P = r.C + t *r.u;
		Vector AB = B - A;
		Vector AC = C - A;
		Vector AP = P - A;
		
		double num = dot(AB,AB) * dot(AC,AC) - dot(AB,AC) * dot(AB,AC);
		

		double beta  = (dot(AC,AC) * dot(AP,AB) - dot(AC,AB) * dot (AP,AC)) / num ;
		//cout << "beta" << endl;
		if (beta <0){return false;}
		if (beta >1){return false;}
		double gamma = (dot(AB,AB) * dot(AP,AC) - dot(AB,AC) * dot (AP,AB)) / num ; 
		if (gamma <0){return false;}
		if (gamma >1){return false;}
		if (beta + gamma > 0.999){return false;}
		//cout << "true" << endl;

		return true;
	}
	Vector A; Vector B; Vector C;
};


class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class Bbox{
public:
	Bbox(){vmin = Vector(1E9,1E9,1E9);vmax = Vector(-1E9,-1E9,-1E9);}
	
	bool inter(const Ray& r){
		double t1x = (vmin.x - r.C.x) / r.u.x;
		double t2x = (vmax.x - r.C.x) / r.u.x;
		double txmin = min(t1x, t2x);
		double txmax = max(t1x, t2x);
		double t1y = (vmin.y - r.C.y) / r.u.y;
		double t2y = (vmax.y - r.C.y) / r.u.y;
		double tymin = min(t1y, t2y);
		double tymax = max(t1y, t2y);
		double t1z = (vmin.z - r.C.z) / r.u.z;
		double t2z = (vmax.z - r.C.z) / r.u.z;
		double tzmin = min(t1z, t2z);
		double tzmax = max(t1z, t2z);
		
		
		return ( min(min(txmax,tymax),tzmax) >= max(max(txmin,tymin),tzmin) );
	}
	Vector vmin, vmax;
};

class Geometry :public Object {
public:
  ~Geometry() {}
	Geometry( char* filename,  Vector trans, double scale, Vector albedo, bool mirror, bool transparent, double nsphere, double emissivite){
		this->albedo = albedo;
		this->mirror = mirror;
		this->transparent = transparent;
		this->nsphere = nsphere;
		this->emissivite = emissivite;
		readOBJ(filename);
		for (int i=0; i<indices.size(); i++){
			vertices[i] = scale * Vector(vertices[i].x,vertices[i].z,vertices[i].y) + trans;}
			
	
		Vector vmin = bbox.vmin;
		Vector vmax = bbox.vmax;
		for (int i=0; i<indices.size(); i++){
				vmin.x = min(vertices[i].x,vmin.x);
				vmin.y = min(vertices[i].y,vmin.y);
				vmin.z = min(vertices[i].z,vmin.z);
				
				vmax.x = max(vertices[i].x,vmax.x);
				vmax.y = max(vertices[i].y,vmax.y);
				vmax.z = max(vertices[i].z,vmax.z);
				}
			
		bbox.vmin = vmin;
		bbox.vmax = vmax;		
	}
	void readOBJ(const char* obj) {
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec.x, &vec.y, &vec.z, &col.x, &col.y, &col.z) == 6) {
					col.x = std::min(1., std::max(0., col.x));
					col.y = std::min(1., std::max(0., col.y));
					col.z = std::min(1., std::max(0., col.z));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);  
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec.x, &vec.y, &vec.z); 
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec.x, &vec.y);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}
	
	
	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t){
		if (! bbox.inter(r)){ return false;}
		
		double smallestt = 1E15;
		bool has_inter = false;

		for (int i=0; i<indices.size(); i++){
			Triangle trilocal(vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk], Vector(1,1,1), false, false,12,1);
			Vector Plocal, Nlocal;
			double tlocal;
			bool inter = trilocal.intersect(r, Plocal, Nlocal,tlocal);
			if (inter && tlocal < smallestt){smallestt = tlocal; t = tlocal;P = Plocal; N=Nlocal;}
			if (inter){has_inter=true;}
		}
		return has_inter;
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	Bbox bbox;
	
};

//class BVHNode {
//public:
//	BVHNode() {};
//	
//	Bbox b;
//	int debut, fin;
//};

class Source {
public:
	Source(const Vector& V,double rayon, double I): V(V), rayon(rayon), I(I) {};
	Vector V;
	double rayon;
	double I;
	
};




std::default_random_engine e;


std::uniform_real_distribution<double> u(0.,1.);
int nbrebondMAX = 2;

class Scene {
public: 
	Scene() {};
	
	void addSphere( Object* s) {spheres.push_back(s);}
		
	bool intersection(const Ray& r, Vector &P, Vector &N, double &t, int &object){
		double smallestt = 1E15;
		bool has_inter = false;
		for (int i=0; i<spheres.size(); i++){
			Vector Plocal, Nlocal;
			double tlocal;
			bool inter = spheres[i]->intersect(r, Plocal, Nlocal,tlocal);
			if (inter && tlocal < smallestt){smallestt = tlocal; object = i; t = tlocal;P = Plocal; N=Nlocal;}
			if (inter){has_inter =true;}
		}
		return has_inter;
	}
	
	
	Vector random_cos(const Vector& N){
		double r1 = u(e);
		double r2 = u(e);
		double wx = cos(2*M_PI*r1)*sqrt(1-r2);
		double wy = sin(2*M_PI*r1)*sqrt(1-r2);
		double wz = sqrt(r2);
		Vector v(wx,wy,wz);
		
		Vector T1;
		Vector T2;
		if (N.x <= N.y && N.x <= N.z){T1 = Vector(0,-N.z,N.y);}
		else if (N.y <= N.x && N.y <= N.z){T1 = Vector(-N.z,0,N.x);}
		else {T1 = Vector(-N.y,N.x,0);}
		T1.normalize();
		T2 = cross(N,T1);
		return  v.x * T1 + v.y *T2 + v.z * N;
		
	}
	
	Vector eclairageDirect(Source sourceLumineuse, Vector Pinc, Vector Normale, Vector albedo){
		Vector O = sourceLumineuse.V;
		Vector OX = Pinc - O;
		Vector OXnormal = OX;
		OXnormal.normalize();
		Vector OXprime = sourceLumineuse.rayon * random_cos(OXnormal);
		Vector Xprime = O + OXprime;
		Vector XXprime = Xprime - Pinc;
		bool VXXprime = dot(Normale,OXprime) < 0;
		Vector omeagaI = XXprime;
		double XX2 = XXprime.norm2();
		omeagaI.normalize();Normale.normalize();XXprime.normalize();OXprime.normalize();
		Vector direct;
		
		// On vérifie qu'il n'y a pas d'intersection avant la lumière.
		Ray rayon_lumiere(Pinc + 0.001*Normale,XXprime);
		Vector Pprime, Nprime;int objectprime;double tprime;
		bool ombre = intersection(rayon_lumiere, Pprime,Nprime, tprime,objectprime);
		if (ombre && tprime*tprime < XX2){direct = Vector(0.,0.,0.);} 
		else{direct = sourceLumineuse.I/(4*M_PI_2)  * (dot(omeagaI,Normale)*dot(XXprime,Normale)*(VXXprime?1.:0.)/(XX2*dot(OXnormal,OXprime))) * albedo;}
		return direct;		
		}
	
	
	Vector getColor(const Ray& ray, int nbrebonds, Source L){
		if (nbrebonds < 0){return Vector(0,0,0);}
		//Source L(Vector(-10, 20, 40),3.,2000000000);
		Vector rayColor;
		Vector P, N;double t;int object;
		Sphere Spherelum(L.V, L.rayon, Vector(1,1,1), false, false,1,1);
		if (Spherelum.intersect(ray,P,N,t)){
			if (nbrebonds == nbrebondMAX){return L.I *  Vector(1,1,1);}
			else {return Vector(0,0,0);}
			}
		if (intersection(ray,P,N,t,object)){
			// MIRROIR
			if (spheres[object]->mirror ){
					Vector R = ray.u - 2 * dot(ray.u, N)*N;
					Ray rayR(P + 0.001 * N, R);
					rayColor = getColor(rayR,nbrebonds-1,L);
					}
			
			
			//TRANSPARENT
			
			else if (spheres[object]->transparent){
				double nair = 1.;
				double rapportN = (nair/spheres[object]->nsphere);
				double incidence = dot(ray.u,N);
				Vector Nused = N;
				
				if (incidence > 0){rapportN = 1./rapportN; Nused = - N;}
				
				double racine = 1 -  rapportN * rapportN * (1-dot(ray.u, Nused)*dot(ray.u, Nused)); 
				if (racine > 0){
					Vector R =  rapportN*ray.u - (rapportN * dot(ray.u,Nused) + sqrt(racine))*Nused;
					Ray rayR(P - 0.001 * Nused, R);
					rayColor = getColor(rayR,nbrebonds-1,L);
				}
			}
		
			else {
				//Eclairage direct
				rayColor = rayColor +  eclairageDirect( L,  P, N, spheres[object]->albedo);
						
				//Eclairage indirect
				Vector reflechi = random_cos(N);
				Ray ray_reflechi(P+0.001*N, reflechi);
				rayColor = rayColor + spheres[object]->albedo * getColor(ray_reflechi, nbrebonds-1,L);
			}
		}
		else {rayColor = Vector(0.,0.,0.);}

		return rayColor;
		
	}
	
	//std::vector<Sphere> spheres;
	std::vector<Object *> spheres;
};



double mu = 0.;
double sigma = .2;


std::normal_distribution<double> n(mu, sigma);


int main() {
	cout << "HELLO"<< endl;
	

    int W = 128;
    int H = 128;
    int Nray = 4;
    
    double fov = 60 * M_PI /180.0;
    double tanoffov = tan(fov*0.5);
    
    Sphere s(Vector(8,-9,0) , 5, Vector(1,1,1), true, false,12,1);
    Sphere sbis(Vector(-8,-9,0) , 5, Vector(1,1,1), false, false,2.5,1);
    Sphere stris(Vector(0,-9,20) , 2, Vector(1,0,1), false, true,1.4,1.4);
    
    Triangle triangle(Vector(8.2,-9.4,0.3), Vector(-18,-9,0), Vector(0,-9,20) , Vector(1,0,1), false, false,12,1);
    
    Source Slum(Vector(0, 20, 30),7,2000000000);
    
    Geometry g("BeautifulGirl.obj", Vector(1,0,1), 15, Vector(0,-10,0), false, false,12,1);
    
    Sphere plafond(Vector(0,1000,0) , 940, Vector(.90,.20,.10),false, false,1,1);
    Sphere s3(Vector(0,0,1000) , 940, Vector(.10,.90,.90),false, false,1,1);
    Sphere s4(Vector(0,0,-1000) , 940, Vector(.10,.80,.230),false, false,1,1);
	Sphere sol(Vector(0,-1000,0) , 980, Vector(.90,.20,.10),false, false,1,1);
    Sphere s6(Vector(1000,0,0) , 940, Vector(.56,.23,.56),false, false,1,1);
	Sphere s7(Vector(-1000,0,0) , 960, Vector(.90,.99,.10),false, false,1,1);
    Scene scene;
    scene.addSphere(&s);
    //scene.addSphere(&sbis);
    //scene.addSphere(&stris);
    
    scene.addSphere(&g);
    
    //scene.addSphere(&triangle);
    scene.addSphere(&plafond);
    scene.addSphere(&s3);
    scene.addSphere(&s4);
    scene.addSphere(&sol);
    scene.addSphere(&s6);
	scene.addSphere(&s7);
  
    
    Vector C(0,0,55);
    
    std::vector<unsigned char> image(W*H * 3, 0);
    
    
    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
		cout << i << endl;
        for (int j = 0; j < W; j++) {
			

			Vector pixelColor(0,0,0);
			for (int nr = 0; nr < Nray; nr++){
					// Rayon aléatoire généré dans le pixel 
					double r1 = u(e);
					double r2 = u(e);
					double offsetx = cos(2*M_PI*r1)*log(1-log(r2));
					double offsety = sin(2*M_PI*r1)*sqrt(1-log(r2));
					
					// Rotation de la camera
					double theta = 0;//M_PI/8;
					Vector direction(0,sin(theta),cos(theta));
					Vector up(0,- sin(theta - M_PI/2),cos(theta - M_PI/2));
					Vector right = cross(direction,up);
					
					Vector u(j - W/2 + offsetx, H/2-i + offsety, - W /(2 * tanoffov));
					Vector v = (j - W/2 + offsetx)*right +  (H/2-i + offsety)*up + (- W /(2 * tanoffov))*direction;
					u.normalize();		
					v.normalize();
					//Ray ray(C,u);
					
					// Le rayon prêt à partir
					Ray ray(C,v);
				
				//Le rayon part 
				pixelColor = pixelColor + scene.getColor(ray, nbrebondMAX,Slum);}
				
			// Moyenne des rayons partis
			pixelColor = pixelColor / Nray;


			//Mise à jour des pixels
			image[(i*W + j) * 3 + 0] = std::min(255., pow(pixelColor.x,0.45));
            image[(i*W + j) * 3 + 1] = std::min(255., pow(pixelColor.y,0.45));
            image[(i*W + j) * 3 + 2] = std::min(255., pow(pixelColor.z,0.45));
        }
	 
    }
    cout << "BYE BYE";
    stbi_write_png("bboxsmall.png", W, H, 3, &image[0], 0);

    return 0;
}

/*
Surface transparente
Transmission Loi de descartes
*/
