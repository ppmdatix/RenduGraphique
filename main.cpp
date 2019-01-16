#include <iostream>

using namespace std;

#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "algorithm"


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

// A ray is an origin and a direction
class Ray {
public:
	Ray(const Vector& C , const Vector& u): C(C), u(u) {};
	Vector C;
	Vector u;
}	;
class Sphere {
public:
    Sphere(const Vector& O, double R, Vector albedo, bool mirror, bool transparent): O(O), R(R), albedo(albedo), mirror(mirror), transparent(transparent) {};
    
    bool intersect(const Ray& r, Vector &P, Vector &N, double &t){
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
		return (delta>=0);}
	
    Vector O;
    double R;
    Vector albedo;
    bool mirror;
    bool transparent;
};


class Source {
public:
	Source(const Vector& V, double I): V(V), I(I) {};
	Vector V;
	double I;
	
};


class Scene {
public: 
	Scene() {};
	
	void addSphere(const Sphere& s) {spheres.push_back(s);}
		
	bool intersection(const Ray& r, Vector &P, Vector &N, double &t, int &object){
		double smallestt = 1E15;
		bool has_inter = false;
		
		for (int i=0; i<spheres.size(); i++){
			Vector Plocal, Nlocal;
			double tlocal;
			bool inter = spheres[i].intersect(r, Plocal, Nlocal,tlocal);
			if (inter && t < smallestt){smallestt = tlocal; object = i; t = tlocal;P = Plocal; N=Nlocal;}
			if (inter){has_inter =true;}
			}
		return has_inter;
		}
		
		
	
	
	Vector getColor(const Ray& ray, int nbrebonds){
		Source L(Vector(-10, 20, 40),2000000000);
		Vector rayColor;
		Vector P, N;double t;int object;
		if (intersection(ray,P,N,t,object)){
			if (spheres[object].mirror && nbrebonds > 0){
					Vector R = ray.u - 2 * dot(ray.u, N)*N;
					Ray rayR(P + 0.001 * N, R);
					rayColor = getColor(rayR,nbrebonds-1);
					}
			else {
				Vector PL = L.V - P;
				double distanceLum2 = PL.norm2();
				PL.normalize();
				Ray rayon_lumiere(P + 0.001*N,PL);
				Vector Pprime, Nprime;int objectprime;double tprime;
				bool ombre = intersection(rayon_lumiere, Pprime,Nprime, tprime,objectprime);
				if (ombre && tprime*tprime < distanceLum2){rayColor = Vector(0.,0.,0.);} 
				else {
					Sphere sphere2 = spheres[object];
					rayColor.x = L.I * (sphere2.albedo.x / M_PI) * max(0.,dot(N,PL)) / (4 * M_PI * distanceLum2);
					rayColor.y = L.I * (sphere2.albedo.y / M_PI) * max(0.,dot(N,PL)) / (4 * M_PI * distanceLum2);
					rayColor.z = L.I * (sphere2.albedo.z / M_PI) * max(0.,dot(N,PL)) / (4 * M_PI * distanceLum2);
					 }
				}
			
			}
		else {rayColor = Vector(0.,0.,0.);}

		return rayColor;
		
		}
	
	std::vector<Sphere> spheres;
};






int main() {
    int W = 512;
    int H = 512;
    
    double fov = 60 * M_PI /180.0;
    double tanoffov = tan(fov*0.5);
    
    Sphere s(Vector(0,0,0) , 10, Vector(1,1,1), false, false);
    Sphere s2(Vector(0,1000,0) , 940, Vector(.90,.20,.10),false, false);
    Sphere s3(Vector(0,0,1000) , 940, Vector(.10,.90,.90),false, false);
    Sphere s4(Vector(0,0,-1000) , 940, Vector(.10,.80,.230),false, false);
	Sphere s5(Vector(0,-1000,0) , 1010, Vector(.90,.20,.10),false, false);
    Sphere s6(Vector(1000,0,0) , 940, Vector(.56,.23,.56),false, false);
	Sphere s7(Vector(-1000,0,0) , 960, Vector(.90,.99,.10),false, false);
    Scene scene;
    scene.addSphere(s);
    scene.addSphere(s2);
    scene.addSphere(s3);
    scene.addSphere(s4);
    scene.addSphere(s5);
    scene.addSphere(s6);
	scene.addSphere(s7);
    
    
    //Source L(Vector(-10, 20, 40),500000000000);
    
    Vector C(0,0,55);
    
    
    
    std::vector<unsigned char> image(W*H * 3, 0);
    
    
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
			
			Vector u(j - W/2, H/2-i, - W /(2 * tanoffov));
			u.normalize();		
			Ray ray(C,u);
			Vector pixelColor;
			pixelColor = scene.getColor(ray, 5);

			image[(i*W + j) * 3 + 0] = std::min(255., pow(pixelColor.x,0.45));
            image[(i*W + j) * 3 + 1] = std::min(255., pow(pixelColor.y,0.45));
            image[(i*W + j) * 3 + 2] = std::min(255., pow(pixelColor.z,0.45));
            
        }
	 
    }
    stbi_write_png("image3.png", W, H, 3, &image[0], 0);

    return 0;
}

/*
Surface transparente
Transmission Loi de descartes
*/
