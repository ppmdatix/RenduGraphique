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
	return Vector(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
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
    Sphere(const Vector& O, double R, Vector albedo, bool mirror, bool transparent, double nsphere, double emissivite): O(O), R(R), albedo(albedo), mirror(mirror), transparent(transparent), nsphere(nsphere), emissivite(emissivite){};
    
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
    double emissivite;
};


class Source {
public:
	Source(const Vector& V,double rayon, double I): V(V), rayon(rayon), I(I) {};
	Vector V;
	double rayon;
	double I;
	
};


std::default_random_engine e;


std::uniform_real_distribution<double> u(0.,1.);
int nbrebondMAX = 5;

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
			if (spheres[object].mirror ){
					Vector R = ray.u - 2 * dot(ray.u, N)*N;
					Ray rayR(P + 0.001 * N, R);
					rayColor = getColor(rayR,nbrebonds-1,L);
					}
			
			
			//TRANSPARENT
			
			else if (spheres[object].transparent){
				double nair = 1.;
				double rapportN = (nair/spheres[object].nsphere);
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
				rayColor = rayColor +  eclairageDirect( L,  P, N, spheres[object].albedo);
						
				//Eclairage indirect
				Vector reflechi = random_cos(N);
				Ray ray_reflechi(P+0.001*N, reflechi);
				rayColor = rayColor + spheres[object].albedo * getColor(ray_reflechi, nbrebonds-1,L);
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




//#pragma omp parallel for



int main() {

    int W = 512;
    int H = 512;
    int Nray = 100;
    
    double fov = 60 * M_PI /180.0;
    double tanoffov = tan(fov*0.5);
    
    Sphere s(Vector(8,-9,0) , 5, Vector(1,1,1), true, false,12,1);
    Sphere sbis(Vector(-8,-9,0) , 5, Vector(1,1,1), false, false,2.5,1);
    Sphere stris(Vector(0,-9,20) , 2, Vector(1,0,1), false, true,1.4,1.4);
    Source Slum(Vector(0, 20, 30),7,2000000000);
    
    
    Sphere plafond(Vector(0,1000,0) , 940, Vector(.90,.20,.10),false, false,1,1);
    Sphere s3(Vector(0,0,1000) , 940, Vector(.10,.90,.90),false, false,1,1);
    Sphere s4(Vector(0,0,-1000) , 940, Vector(.10,.80,.230),false, false,1,1);
	Sphere sol(Vector(0,-1000,0) , 980, Vector(.90,.20,.10),false, false,1,1);
    Sphere s6(Vector(1000,0,0) , 940, Vector(.56,.23,.56),false, false,1,1);
	Sphere s7(Vector(-1000,0,0) , 960, Vector(.90,.99,.10),false, false,1,1);
    Scene scene;
    scene.addSphere(s);
    scene.addSphere(sbis);
    scene.addSphere(stris);
    scene.addSphere(plafond);
    scene.addSphere(s3);
    scene.addSphere(s4);
    scene.addSphere(sol);
    scene.addSphere(s6);
	scene.addSphere(s7);
  
    
    Vector C(0,0,55);
    
    std::vector<unsigned char> image(W*H * 3, 0);
    
    
    

    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
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
    stbi_write_png("image30_01.png", W, H, 3, &image[0], 0);

    return 0;
}

/*
Surface transparente
Transmission Loi de descartes
*/
