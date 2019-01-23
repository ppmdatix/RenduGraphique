#include <iostream>

using namespace std;

#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>
#include <random>

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
	return Vector(a.y*b.z-a.z*b.y, a.y*b.z-a.z*b.y,a.y*b.z-a.y*b.x);
}

// A ray is an origin and a direction
class Ray {
public:
	Ray(const Vector& C , const Vector& u): C(C), u(u) {};
	Vector C;
	Vector u;
}	;
class Sphere {
public:
    Sphere(const Vector& O, double R, Vector albedo, bool mirror, bool transparent, double nsphere): O(O), R(R), albedo(albedo), mirror(mirror), transparent(transparent), nsphere(nsphere){};
    
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
		return (delta>=0);
	}
	
    Vector O;
    double R;
    Vector albedo;
    bool mirror;
    bool transparent;
    double nsphere;
};


class Source {
public:
	Source(const Vector& V, double I): V(V), I(I) {};
	Vector V;
	double I;
	
};


std::default_random_engine e;


std::uniform_real_distribution<double> u(0.,1.);

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
			else if (spheres[object].transparent && nbrebonds > 0){
				double rapportN = (1/spheres[object].nsphere);
				
				double racine =1 -  rapportN * rapportN * (1-dot(ray.u, N)*dot(ray.u, N)); 
				if (racine < 0){
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
				else {
					Vector R =  rapportN*ray.u - (rapportN * dot(ray.u,N) + sqrt(racine))*N;
					Ray rayR(P + 0.001 * N, R);
					rayColor = getColor(rayR,nbrebonds-1);
				
					
				}
					
			}
		
			else {
				
				//Eclairage direct
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
				//Eclairage indirect
				
				Vector reflechi = random_cos(N);
				Ray ray_reflechi(P+0.001*N, reflechi);
				rayColor = rayColor + spheres[object].albedo * getColor(ray_reflechi, nbrebonds-1);

				
				
			}
		}
		else {rayColor = Vector(0.,0.,0.);}

		return rayColor;
		
	}
	
	std::vector<Sphere> spheres;
};




double mu = 0.;
double sigma = .2;

std::normal_distribution<double> n(mu, sigma);

class integrercos10 {
public:
	integrercos10(int N): N(N) {};
	double generate(){
		double somme = 0;
		
		for  (int i=0; i<N; i++){
			double theta = -M_PI/2 + u(e) * M_PI;
			somme = somme + pow(cos(theta),10.) * M_PI;
		}
		somme = somme / N;
		cout << somme;
		
		return somme;
		}
		
	double generateNormal(){
		double somme = 0;
		
		for  (int i=0; i<N; i++){
			double theta = n(e);
			if (abs(theta)<(M_PI/2)){
				double ptheta = exp(-pow(theta-mu,2)/(2*sigma*sigma))/(sigma*sqrt(2+M_PI));
				somme = somme + pow(cos(theta),10.) / ptheta;
				}
			else {cout<< "missed";}
			//somme = somme + pow(cos(theta),10.) * M_PI;
		}
		somme = somme / N;
		cout << somme;
		
		return somme;
		}
	double N;
	
};

class integrercoscoscos {
public:
	integrercoscoscos(int N): N(N) {};

		
	double generateNormal(){
		double somme = 0;
		
		for  (int i=0; i<N; i++){
			double theta1 = n(e);
			double theta2 = n(e);
			double theta3 = n(e);
			if (abs(theta1)<(M_PI/2) and abs(theta2)<(M_PI/2) and abs(theta3)<(M_PI/2)){
				double ptheta1 = exp(-pow(theta1-mu,2)/(2*sigma*sigma))/(sigma*sqrt(2+M_PI));
				double ptheta2 = exp(-pow(theta2-mu,2)/(2*sigma*sigma))/(sigma*sqrt(2+M_PI));
				double ptheta3 = exp(-pow(theta3-mu,2)/(2*sigma*sigma))/(sigma*sqrt(2+M_PI));
				
				somme = somme + pow(cos(theta1*theta2),30.)*pow(cos(theta3*theta2),30.) / (ptheta1*ptheta2*ptheta3);
				}
			else {cout<< "missed";}
			//somme = somme + pow(cos(theta),10.) * M_PI;
		}
		somme = somme / N;
		cout << somme;
		
		return somme;
		}
	double N;
	
};






int main() {
	
	integrercoscoscos inte(10000);
	inte.generateNormal();
    int W = 512;
    int H = 512;
    int Nray = 5;
    
    double fov = 60 * M_PI /180.0;
    double tanoffov = tan(fov*0.5);
    
    Sphere s(Vector(8,0,0) , 5, Vector(1,1,1), true, true,12);
    Sphere sbis(Vector(-8,0,0) , 5, Vector(1,0,1), false, false,1);
    
    
    Sphere s2(Vector(0,1000,0) , 940, Vector(.90,.20,.10),false, false,1);
    Sphere s3(Vector(0,0,1000) , 940, Vector(.10,.90,.90),false, false,1);
    Sphere s4(Vector(0,0,-1000) , 940, Vector(.10,.80,.230),false, false,1);
	Sphere s5(Vector(0,-1000,0) , 960, Vector(.90,.20,.10),false, false,1);
    Sphere s6(Vector(1000,0,0) , 940, Vector(.56,.23,.56),false, false,1);
	Sphere s7(Vector(-1000,0,0) , 960, Vector(.90,.99,.10),false, false,1);
    Scene scene;
    scene.addSphere(s);
    scene.addSphere(sbis);
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
			Vector pixelColor(0,0,0);
			for (int nr = 0; nr < Nray; nr++){pixelColor = pixelColor + scene.getColor(ray, 3);}
			pixelColor = pixelColor / Nray;

			image[(i*W + j) * 3 + 0] = std::min(255., pow(pixelColor.x,0.45));
            image[(i*W + j) * 3 + 1] = std::min(255., pow(pixelColor.y,0.45));
            image[(i*W + j) * 3 + 2] = std::min(255., pow(pixelColor.z,0.45));
        }
	 
    }
    stbi_write_png("imageNeige.png", W, H, 3, &image[0], 0);

    return 0;
}

/*
Surface transparente
Transmission Loi de descartes
*/
