#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include "limits"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 256
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

struct Intersection
{
	vec4 position;
	float distance;
	int triangleIndex;
};

float focalLength = 256;
vec4 cameraPos(0,0,-3,1.0);
float yaw = 0;


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
mat4 setRotationMat(mat4 R,float yaw);
bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection );
//void LoadTestModel( std::vector<Triangle>&triangles );

vector<Triangle> triangles;


//void ray_direction(int x, int y);

int main( int argc, char* argv[] )
{
	
	screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
	LoadTestModel(triangles);

	while ( Update())
		{
			Draw(screen);
			SDL_Renderframe(screen);
		}

	SDL_SaveImage( screen, "screenshot.bmp" );

	KillSDL(screen);
	return 0;
}

/*Place your drawing here*/
void Draw(screen* screen)
{

	memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

	Intersection closestIntersection;


	for(int y = 0; y < SCREEN_HEIGHT; y++) {

		for(int x = 0; x < SCREEN_WIDTH; x++) {
			vec4 d = vec4(x - (SCREEN_WIDTH / 2), y - (SCREEN_HEIGHT / 2), focalLength, 1.0);
			bool b = ClosestIntersection(cameraPos,d, triangles, closestIntersection);
			if(b) {
				int i = closestIntersection.triangleIndex;
				vec3 colour = triangles[i].color;
				PutPixelSDL(screen, x, y, colour);
			} else {
				PutPixelSDL(screen, x, y, vec3(0.0,0.0,0.0));
			}
		}
	}
}
mat4 setRotationMat(mat4 R,float yaw){
	//First Column
	// R = mat4(right, up, forward, -from)
	// R[3][3] = 1
	R[0][0] = cos(yaw);
	R[1][0] = 0;
	R[2][0] = -sin(yaw);
	R[3][0] = 0;
	//Second Column
	R[0][1] = 0;
	R[1][1] = 1;
	R[2][1] = 0;
	R[3][1] = 0;
	//Third Column
	R[0][2] = sin(yaw);
	R[1][2] = 0;
	R[2][2] = cos(yaw);
	R[3][2] = 0;

	//Fourth column
	for(int i=0; i< 4;i++){
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
	static int t = SDL_GetTicks();
	/* Compute frame time */
	int t2 = SDL_GetTicks();
	// float dt = float(t2-t);
	t = t2;
	float change =0.01f;

	// mat4 R;

	float yaw = 0;

	// R = setRotationMat(R,yaw);
	vec4 forward;
	vec4 left;
	vec4 right;
	vec4 down;


	
	// vec4 down(    R[1][0], R[1][1], R[1][2], 1 );
	





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
					case SDLK_UP:
					{
						/* Move camera forward */
						
						// vec4 forward( R[2][0], R[2][1], R[2][2], 1 );
						// cameraPos = forward * vec4(cameraPos[0],cameraPos[1],cameraPos[2],cameraPos[3]);
						cameraPos[2] = cameraPos[2] + change;
						yaw += 0.6f;
						break;
					}
					case SDLK_DOWN:
					{
						/* Move camera backwards */
						// setRotationMat(R,yaw);
						// vec4 down(    R[1][0], R[1][1], R[1][2], 1 );
						// cameraPos = down * vec4(cameraPos[0],cameraPos[1],cameraPos[2],cameraPos[3]);
						cameraPos[2] = cameraPos[2] - change;
						yaw -= 0.6f;
						break;
					}
					case SDLK_LEFT:
					{
						/* Move camera left */
						// setRotationMat(R,yaw);
						// vec4 left(   R[0][0], R[0][1], R[0][2], 1 );
						// cameraPos = left * vec4(cameraPos[0],cameraPos[1],cameraPos[2],cameraPos[3]);
						cameraPos[0] = cameraPos[0] - change;
						yaw -= 0.6f;
						break;
					}
					case SDLK_RIGHT:
					{
						/* Move camera right */
						// setRotationMat(R,yaw);
						// vec4 right(   R[0][0], R[0][1], R[0][2], 1 );
						// cameraPos = right * vec4(cameraPos[0],cameraPos[1],cameraPos[2],cameraPos[3]);
						cameraPos[0] = cameraPos[0] + change;
						yaw += 0.6f;
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



bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection) {

	bool intersects_triangle = false;
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

	 mat3 A2(-d,b,e2); // u 
	 float u = glm::determinant(A2) / glm::determinant(A);

	 mat3 A3(-d, e1, b); // v
	 float v = glm::determinant(A3) / glm::determinant(A);
	 //vec3 dist = t * d;
	 //float mag = sqrt(dist.length());

	// Re--jig Code so only do computation if its close


	 if(t >= 0 && u >= 0 && v >= 0 && (u+v) <= 1) {
		intersects_triangle = true;

		//updating closest intersection
		if(t < closestIntersection.distance) {
			closestIntersection.position = start + (t*dir);
			closestIntersection.triangleIndex = i;
			closestIntersection.distance = t;
		}

	 }
	}
	return intersects_triangle;
}


