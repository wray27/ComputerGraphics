#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include "limits"
#include "ObjLoader.h"
#include <omp.h>



using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;
#define PI 3.14159265358979323546
#define SCREEN_WIDTH 512
#define SCREEN_HEIGHT 512
#define FULLSCREEN_MODE false

struct Intersection
{
	vec4 position;
	float distance;
	int triangleIndex;
	bool is_sphere;
};

struct Light
{
	vec4 pos;
	vec3 col;
};

struct Sphere
{
	vec4 center;
	float radius;
	vec3 colour;
	bool is_reflective;
};

float focalLength = 512;
vec4 cameraPos(0,0,-3,1.0);
float yaw = 0;
//vec4 lightPos = vec4( 0, -0.5, -0.7, 1.0 );




/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
mat4 setRotationMat(mat4 R,float yaw);
bool Update();
vec3 DirectLight( const Intersection& i );
vec3 IndirectLight( const Intersection& i);
vec3 light( const Intersection& i);
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection );
bool raySphereIntersection(vec4 start, vec4 dir, Sphere &s, Intersection &closestInt);
void setupLighting(vector<Light> &lights);
void setupSphere(Sphere &s);
bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1);

//void LoadTestModel( std::vector<Triangle>&triangles );

vector<Triangle> triangles;
vector<Light> lights;
Sphere sphere;

//initialise a light here


/*struct Triangle
   {
       vec4 v0;
       vec4 v1;
       vec4 v2;
       vec4 normal;
       vec3 color;
       //bool is_reflective
};*/



int main( int argc, char* argv[] )
{
	setupLighting(lights);
	screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
	//load in sphere and triangles
	setupSphere(sphere);
	LoadTestModel(triangles);
	//Object Loading
	//int objectCount= 0;
	// if(loadOBJ("./Objects/bunny.obj",triangles)) {
	// 	cout << "Object " << objectCount << " loaded successfully!\n";
	// 	cout << triangles[0].v1.x << "," << triangles[0].v1.y << "," << triangles[0].v1.z << endl;
	// 	objectCount ++;
	// }

	while (Update())
	{


		Draw(screen);
		SDL_Renderframe(screen);
		SDL_SaveImage( screen, "screenshot.bmp" );

	}

	SDL_SaveImage( screen, "screenshot.bmp" );

	KillSDL(screen);
	return 0;
}

void setupLighting(vector<Light> &lights) {
	// function to add more lights if needed
	 lights.resize(1);
	//original light
	 lights[0].pos = vec4( 0, -0.5, -0.7, 1.0 );
	 lights[0].col = 14.f * vec3( 1.0f, 1.0f, 1.0f );
}
void setupSphere(Sphere &s) {
		s.center = vec4(-0.5f,0.7f,-0.7f,1.0f);
		s.radius = 0.3f;
		s.colour = vec3(1.0f,0.27f,0.0f);
		s.is_reflective = 0;


}


/*Place your drawing here*/
void Draw(screen* screen) {

	memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
	Intersection closestIntersection;
	Intersection closestIntersection_two;

  #pragma omp parallel for private(closestIntersection , closestIntersection_two)
	for(int y = 0; y < SCREEN_HEIGHT; y++) {
		for(int x = 0; x < SCREEN_WIDTH; x++) {
			vec3 colour = vec3(0.0f);
			//Anti-Aliasing and Depth of Field
			 for (float i = -0.5; i <= 0.5; i+=0.25f) {
					for (float j = -0.5; j <= 0.5; j +=0.25f) {
				    //vec4 d = vec4(x - (SCREEN_WIDTH / 2), y - (SCREEN_HEIGHT / 2), focalLength, 1.0);
					  vec4 d = vec4(x - (SCREEN_WIDTH / 2) + i*256/20.5f, y - (SCREEN_HEIGHT / 2) + j*256/20.5f, focalLength, 1.0);
					  bool b = ClosestIntersection(cameraPos - vec4(i, j, 0, 0)/20.5f, d, triangles, closestIntersection);
						bool c = raySphereIntersection(cameraPos - vec4(i, j, 0, 0)/20.5f, d,sphere,closestIntersection);
						if(b || c) {
							//int i = closestIntersection.triangleIndex;
							//colour = triangles[i].color;
							colour += IndirectLight(closestIntersection);
						}

					//Reflection
						if(triangles[closestIntersection.triangleIndex].is_reflective) {
							colour = vec3(0.0f,0.0f,0.0f);
							vec4 dir = closestIntersection.position - cameraPos;
							vec4 N = triangles[closestIntersection.triangleIndex].normal;
							//vec4 N = normalize(closestIntersection.position - sphere.center);
							vec4 start = closestIntersection.position;
							vec3 reflected_v3 = glm::reflect(vec3(dir), vec3(N));
							vec4 reflected_ray = vec4(reflected_v3, 1.0f);
							if(ClosestIntersection(start + 1e-3f*reflected_ray - vec4(i, j, 0, 0)/20.5f, reflected_ray, triangles, closestIntersection_two)) {
								colour += IndirectLight(closestIntersection_two);
							}
						}
				  }
			  }
			PutPixelSDL(screen, x, y, colour /25.0f);

		}
		SDL_Renderframe(screen);
	}
}


mat4 setRotationMat(mat4 R, float yaw){
	//First Column
	// R = mat4(right, up, forward, -from)
	// R[3][3] = 1

	float A = 0;
	float B = yaw;
	float C = 0;
	// Compute rotation matrix
	R[0][0] =  cos(C)*cos(B);
	R[1][0] =  sin(C)*cos(B);
	R[2][0] = -sin(B);

	R[0][1] = -sin(C)*cos(A) + cos(C)*sin(B)*sin(A);
	R[1][1] =  cos(C)*cos(A) + sin(C)*sin(B)*sin(A);
	R[2][1] =  cos(B)*sin(A);

	R[0][2] =  sin(C)*sin(A) + cos(C)*sin(B)*cos(A);
	R[1][2] = -cos(C)*sin(A) + sin(C)*sin(B)*cos(A);
	R[2][2] =  cos(B)*cos(A);


	// R[0][0] = cos(yaw);
	// R[1][0] = 0;
	// R[2][0] = -sin(yaw);
	R[3][0] = 0;

	// //Second Column
	// R[0][1] = 0;
	// R[1][1] = 1;
	// R[2][1] = 0;
	R[3][1] = 0;
	// //Third Column
	// R[0][2] = sin(yaw);
	// R[1][2] = 0;
	// R[2][2] = cos(yaw);
	R[3][2] = 0;

	//Fourth column
	for (int i = 0; i < 3; i++){
		R[i][3] = -cameraPos[i];

	}
	// R[0][3] = 1;
	// R[1][3] = 1;
	// R[2][3] = 1;
	R[3][3] = 1;

	return R;

}

/*Place updates of parameters here*/
bool Update()
{
	//vec4 lightPos = lights.at(0).pos;

	static int t = SDL_GetTicks();
	/* Compute frame time */
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	float change =0.1f;
	cout << "Render time: " << dt << "ms." << endl;
	mat4 R;
	R = setRotationMat(R,yaw);
	vec4 forward( R[2][0], R[2][1], R[2][2], 1 );
	vec4 right(   R[0][0], R[0][1], R[0][2], 1 );
  vec4 down(    R[1][0], R[1][1], R[1][2], 1 );



	// if( keystate[SDLK_w] )
    //  lightPos += forward;

	// if( keystate[SDLK_s] )
	// 	lightPos -= forward;

	// if( keystate[SDLK_a] )
	// 	lightPos -= right;

	// if( keystate[SDLK_d] )
	// 	lightPos += right;

	// if( keystate[SDLK_e] )
	// 	lightPos += down;

	// if( keystate[SDLK_q] )

	// 	lightPos -= down;


	SDL_Event e;
	while(SDL_PollEvent(&e))
	{
		if (e.type == SDL_QUIT)
		{
			return false;
		}else{
			if (e.type == SDL_KEYDOWN)
			{
				int key_code = e.key.keysym.sym;
				switch(key_code)
				{
					case SDLK_w:
					{
						lights.at(0).pos += forward;
						break;
					}
					case SDLK_s:
					{
						lights.at(0).pos -= forward;
						break;
					}
					case SDLK_a:
					{
						lights.at(0).pos += right;
						break;
					}
					case SDLK_d:
					{
						lights.at(0).pos -= right;
						break;
					}
					case SDLK_e:
					{
						lights.at(0).pos -= down;
						break;
					}
					case SDLK_q:
					{
						lights.at(0).pos += down;
						break;
					}
					case SDLK_UP:
					{
						/* Move camera forward */

						// vec4 forward( R[2][0], R[2][1], R[2][2], 1 );
						// cameraPos = forward * vec4(cameraPos[0],cameraPos[1],cameraPos[2],cameraPos[3]);
						cameraPos[2] = cameraPos[2] + change;
						yaw += change;
						break;
					}
					case SDLK_DOWN:
					{
						/* Move camera backwards */
						// setRotationMat(R,yaw);
						// vec4 down(    R[1][0], R[1][1], R[1][2], 1 );
						// cameraPos = down * vec4(cameraPos[0],cameraPos[1],cameraPos[2],cameraPos[3]);
						cameraPos[2] = cameraPos[2] - change;
						yaw -= change;
						break;
					}
					case SDLK_LEFT:
					{
						/* Move camera left */
						// setRotationMat(R,yaw);
						// vec4 left(   R[0][0], R[0][1], R[0][2], 1 );
						// cameraPos = left * vec4(cameraPos[0],cameraPos[1],cameraPos[2],cameraPos[3]);
						cameraPos[0] = cameraPos[0] - change;
						yaw -= change;
						break;
					}
					case SDLK_RIGHT:
					{
						/* Move camera right */
						// setRotationMat(R,yaw);
						// vec4 right(   R[0][0], R[0][1], R[0][2], 1 );
						// cameraPos = right * vec4(cameraPos[0],cameraPos[1],cameraPos[2],cameraPos[3]);
						cameraPos[0] = cameraPos[0] + change;
						yaw += change;
						break;
					}
					case SDLK_ESCAPE:
					{
						/* Move camera quit */
						return false;
					}
						// /* Move camera quit */
						// return false;
				}
			}
		}
		return true;
	}
	return true;



}

// added soft shadows
vec3 light( const Intersection& i){

	vec3 colour = vec3(0.0f,0.0f,0.0f);
	vec4 lightPos1 = lights.at(0).pos;
	vec3 lightColor = lights.at(0).col;

	for(float k = -0.2f; k <= 0.2f; k+= 0.05f) {
		for(float j = -0.2f; j <= 0.2f; j+= 0.05f){
			vec4 lightPos = lightPos1 + vec4(k, 0, j, 0);

			//vec4 lightPos = lights.at(j).pos;

			// calculates whether an area has a shadow
			Intersection shadowI;
			vec4 dir = lightPos - i.position;
			ClosestIntersection(i.position+1e-3f*dir,dir,triangles,shadowI);
			raySphereIntersection(i.position+1e-3f*dir,dir,sphere, shadowI);


			//if the distance is greater than 1 then there is not a shadow
			if (shadowI.distance >= 1) {
				//color vector describes the power P
				//vec3 lightColor = 14.f * vec3( 1.0f, 1.0f, 1.0f );

				//calculatiing the distance from light source - r
				float distance = glm::distance(lightPos,i.position);

				// unit vector describing the direction from the surface point to the light source
				vec4 _r = glm::normalize(lightPos - i.position);
				// vec3 r = glm::normalize(vec3(_r.x,_r.y,_r.z));

				//normal pointing out of the surface as a unit vector
				vec4 _n = i.is_sphere ? normalize(i.position - sphere.center) : glm::normalize(triangles[i.triangleIndex].normal);
				// vec3 n = glm::normalize(vec3(triangles[i.triangleIndex].normal.x,triangles[i.triangleIndex].normal.y,triangles[i.triangleIndex].normal.x));

				// Adds light to the scene with no colour
				colour += (lightColor * max(glm::dot(_r,_n),0.0f)) / (float)(4*PI*pow(distance,2));
				// adds colour tot the scene
				// colour = colour * triangles[i.triangleIndex].color;
			}
		}
	}


	colour /= 60.0f;

	return colour;
}


vec3 DirectLight(const Intersection& i ){
	vec3 ret = light(i) * triangles[i.triangleIndex].color;
	//vec3 ret = i.is_sphere ? light(i) * sphere.colour : light(i) * triangles[i.triangleIndex].color;
	return ret;
}


vec3 IndirectLight(const Intersection& i) {
	vec3 colour;
	vec3 indirectLight = 0.5f*vec3(1.0f,1.0f,1.0f);
	// float distance = glm::distance(lightPos,i.position);
	vec3 illum = light(i);
	colour = i.is_sphere ? (illum + indirectLight) * sphere.colour : (illum + indirectLight) * triangles[i.triangleIndex].color;

	// adds colour tot the scene
	//colour = colour * triangles[i.triangleIndex].color;
	//colour = triangles[i.triangleIndex].color * total;

	return colour;
}






bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection) {

	bool intersectsTriangle = false;
	closestIntersection.distance = std :: numeric_limits<float>::max();
	mat4 R;
	R = setRotationMat(R,yaw);

	for(int i = 0; i < triangles.size(); i++) {

	 Triangle triangle = triangles[i];
	 vec4 v0 = R*triangle.v0;
	 vec4 v1 = R*triangle.v1;
	 vec4 v2 = R*triangle.v2;
	 vec3 e1 = vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
	 vec3 e2 = vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
	 vec3 b = vec3(start.x-v0.x,start.y-v0.y,start.z-v0.z);
	 vec3 d = vec3(dir.x, dir.y,dir.z);

	 mat3 A(-d, e1, e2);

	 // vec3 x = glm::inverse( A ) * b;

	 // float t = x.x;
	 // float u = x.y;
	 // float v = x.z;
	 mat3 A1(b, e1, e2); //t
	 float t = glm::determinant(A1) / glm::determinant(A);

	 if(t > 0.001) {

	  mat3 A2(-d,b,e2); // u
	  float u = glm::determinant(A2) / glm::determinant(A);

	  mat3 A3(-d, e1, b); // v
	  float v = glm::determinant(A3) / glm::determinant(A);
	  //vec3 dist = t * d;
	  //float mag = sqrt(dist.length());

	  // Re--jig Code so only do computation if its close


	  if(u >= 0 && v >= 0 && (u+v) <= 1) {
		 intersectsTriangle = true;

		 //updating closest intersection
		 if(t < closestIntersection.distance) {
			closestIntersection.position = start + (t*dir);
			closestIntersection.triangleIndex = i;
			closestIntersection.distance = t;
			closestIntersection.is_sphere = false;
		 }
	  }
	 }
	}

	return intersectsTriangle;
}

bool raySphereIntersection(vec4 start, vec4 dir, Sphere &sphere, Intersection &closestIntersection) {

	float radius2 = (sphere.radius * sphere.radius);
	float t0, t1; //solutions - if there are any
	bool intersects_sphere = false;
	// vec4 L = sphere.center - start;
	// float tca = dot(L,dir);
	// float d2 = dot(L,L) - (tca * tca);
	// if(d2 > radius2) return false;
	// float thc = sqrt(radius2 - d2);
	// t0 = tca - thc;
	// t1 = tca + thc;
	//
	// if(t0 > t1) std::swap(t0, t1);
	// if(t0 < 0) {
	// 	t0 = t1;
	// 	if(t0 < 0) return false;
	// }
	//
	// if(t0 < closestIntersection.distance) {
	// closestIntersection.distance = t0;
	// closestIntersection.position = start + (t0*dir);
	// return true;
	// }

	vec4 L = start - sphere.center;
	float a = glm::dot(dir,dir);
	float b = 2 * glm::dot(dir,L);
	float c = glm::dot(L,L) - radius2;
	if(!solveQuadratic(a,b,c,t0,t1)) return intersects_sphere;

	if(t0 > t1) std::swap(t0,t1);
	if(t0 < 0) {
		t0 = t1;
		if(t0 < 0) return intersects_sphere;
	}


	if(t0 < closestIntersection.distance) {
	closestIntersection.distance = t0;
	closestIntersection.position = start + (t0*dir);
	closestIntersection.is_sphere = true;
	// return true;
	intersects_sphere = true;
	}

	return intersects_sphere;
}


bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
float discr = b * b - 4 * a * c;
if (discr < 0) return false;
else if (discr == 0) {
x0 = x1 = - 0.5 * b / a;
}
else {
float q = (b > 0) ?
-0.5 * (b + sqrt(discr)) :
-0.5 * (b - sqrt(discr));
x0 = q / a;
x1 = c / q;
}

return true;
}
