#include <iostream>

using namespace std;

#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "algorithm"



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


class Ray {
public:
	Ray(const Vector& C , const Vector& u): C(C), u(u) {};
	Vector C;
	Vector u;
}	;
class Sphere {
public:
    Sphere(const Vector& O, double R, double albedo): O(O), R(R), albedo(albedo) {};
    
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
				else {return false;
				}
			}
		}
		N = P- O;
		N.normalize();
		return (delta>=0);}
	
    Vector O;
    double R;
    double albedo;
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
	
	void addSphere(const Sphere& s) {
		spheres.push_back(s);}
		
		bool intersection(const Ray& r, Vector &P, Vector &N){
			double smallestt = 1E15;
			for (int i=0, i<spheres.size(), i++){
				Vector Plocal, Nlocal;
				double t;
				bool inter = spheres[i].intersect(r, Plocal, Nlocal,t){
					if (inter && t < smallestt){
						smallest = t;
						}
					
					}
				
				}
	 std::vector<Sphere> spheres;
	
};






int main() {
    int W = 512;
    int H = 512;
    
    double fov = 60 * M_PI /180.0;
    double tanoffov = tan(fov*0.5);
    
    Sphere s(Vector(0,0,0) , 10, 1);
    Sphere s2(Vector(40,0,0) , 30, 1);
    
    Source L(Vector(-10, 20, 40),10000000);
    
    Vector C(0,0,55);
    
    
    std::vector<unsigned char> image(W*H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
			
			Vector u(j - W/2, i-H/2, - W /(2 * tanoffov));
			u.normalize();		
			
			Ray ray(C,u);
			Vector pixelColor;
			Vector P, N;
			if(s.intersect(ray,P,N)	){
				Vector PL = L.V - P;
				double distanceLum2 = PL.norm2();
				PL.normalize();
				
				pixelColor = L.I * (s.albedo / M_PI) *std::max(0.,dot(N,PL)) / (4 * M_PI * distanceLum2);
			} else { pixelColor = Vector(0,0,0);}
			image[(i*W + j) * 3 + 0] = std::min(255., pixelColor.x);
            image[(i*W + j) * 3 + 1] = std::min(255., pixelColor.y);
            image[(i*W + j) * 3 + 2] = std::min(255., pixelColor.z);
            

        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}
