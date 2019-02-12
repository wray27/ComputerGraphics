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

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

struct Intersection
{
  vec4 position;
  float distance;
  int triangleIndex;
};

float focalLength = 180;
vec4 cameraPos(0,0.2,-2,1.0);



/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

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
        PutPixelSDL(screen, x, y, vec3(1.0,1.0,1.0));
      }
    }
  }
}

/*Place updates of parameters here*/
bool Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  // float dt = float(t2-t);
  t = t2;
  float change =0.3f;

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
            /* Move camera forward */
            cameraPos=vec4(cameraPos[0],cameraPos[1],cameraPos[2]+change,cameraPos[3]);
            break;
          case SDLK_DOWN:
            /* Move camera backwards */
            cameraPos=vec4(cameraPos[0],cameraPos[1],cameraPos[2]-change,cameraPos[3]);
            break;
        	case SDLK_LEFT:
            /* Move camera left */
            cameraPos=vec4(cameraPos[0]-change,cameraPos[1],cameraPos[2],cameraPos[3]);
            break;
        	case SDLK_RIGHT:
            /* Move camera right */
            cameraPos=vec4(cameraPos[0]+change,cameraPos[1],cameraPos[2],cameraPos[3]);
            break;
        	case SDLK_ESCAPE:
            /* Move camera quit */
            return false;
        }
      }  
    }
  return true;
  }
}



bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection ) {

  bool intersects_triangle = false;
  closestIntersection.distance = std :: numeric_limits<float>::max();

  for(int i = 0; i < triangles.size(); i++) {

   Triangle triangle = triangles[i];
   vec4 v0 = triangle.v0;
   vec4 v1 = triangle.v1;
   vec4 v2 = triangle.v2;
   vec3 e1 = vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
   vec3 e2 = vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
   vec3 b = vec3(start.x-v0.x,start.y-v0.y,start.z-v0.z);
   vec3 d = vec3(dir.x, dir.y,dir.z);

   mat3 A(-d, e1, e2);
   vec3 x = glm::inverse( A ) * b;

   float t = x.x;
   float u = x.y;
   float v = x.z;

   //vec3 dist = t * d;
   //float mag = sqrt(dist.length());


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

mat3 cramersRule(mat3 A) {

  float a = A[0][0];
  float b = A[0][1];
  float c = A[0][2];
  float d = A[1][0];
  float e = A[1][1];
  float f = A[1][2];
  float g = A[2][0];
  float h = A[2][1];
  float i = A[2][2];


  float det = ((a * e * i) + (b * f * g) +(c * d * h)) - ((b * d * i) + (h * f * a) + (g * e * c));
  // A = glm::transpose(glm::adjoint(A));

  return (1/det) * A; 

}

// void adjoint(mat3 A){

//   mat3 

// }
