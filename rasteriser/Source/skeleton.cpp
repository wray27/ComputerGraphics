#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::ivec2;
using glm::vec2;

SDL_Event event;



#define SCREEN_WIDTH 256
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

vector<Triangle> triangles;
vec4 cameraPos( 0, 0, -3.001,1 );
// focal length of the camera
float f = SCREEN_WIDTH;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& v, ivec2& p );
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color );
void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen );

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
    /* Clear buffer */
    memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
   
    vec3 colour(1.0,1.0,1.0);
    vector<vec4> triangleVerts(3);
   

    // ivec2 delta = glm::abs( a - b );
    // int pixels = glm::max( delta.x, delta.y ) + 1;
    
    // vector<ivec2> line( pixels );
    // Interpolate( a, b, line );
    
    for(int i=0; i<triangles.size(); i++)
    {
        triangleVerts[0] = triangles[i].v0;
        triangleVerts[1] = triangles[i].v1;
        triangleVerts[2] = triangles[i].v2;

        DrawPolygonEdges(triangleVerts,screen);
       
        
    }
}
void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen ){
    int V = vertices.size();
    vec3 color( 1, 1, 1 );
    // Transform each vertex from 3D world position to 2D image position:
    vector<ivec2> projectedVertices( V );
    for( int i=0; i<V; ++i )
    {
        VertexShader( vertices[i], projectedVertices[i] );
    }
    // Loop over all vertices and draw the edge from it to the next vertex:
    for( int i=0; i<V; ++i )
    {
        int j = (i+1)%V; // The next vertex
        
        DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
    } 

}

void VertexShader( const vec4& v, ivec2& p ){
    vec4 temp;
    temp.x = v.x - cameraPos.x;
    temp.y = v.y - cameraPos.y;
    temp.z = v.z - cameraPos.z;



    p.x = (f * (temp.x/temp.z)) + (SCREEN_WIDTH /2);
    p.y = (f * (temp.y/temp.z)) + (SCREEN_HEIGHT /2);


}

void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result )
{
    int N = result.size();
    float temp = float(max(N-1,1));
    
    vec2 step;

    step.x = (b.x - a.x) / temp;
    step.y = (b.y - a.y) / temp;

    vec2 current( a );
    for( int i=0; i<N; ++i )
    {
       result[i].x = round(current.x);
       result[i].y = round(current.y);
       current.x += step.x;
       current.y += step.y;
    }
}
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color ){
    ivec2 delta;
    delta.x = glm::abs( a.x - b.x );
    delta.y = glm::abs( a.y - b.y );
    int pixels = glm::max( delta.x, delta.y ) + 1;
    
    vector<ivec2> line( pixels );
    Interpolate( a, b, line );
    for(int i =0; i < pixels;i ++){
        PutPixelSDL( screen, line[i].x, line[i].y, color );

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

  SDL_Event e;
  while(SDL_PollEvent(&e))
    {
      if (e.type == SDL_QUIT)
	{
	  return false;
	}
      else
	if (e.type == SDL_KEYDOWN)
	  {
	    int key_code = e.key.keysym.sym;
	    switch(key_code)
	      {
	      case SDLK_UP:
		/* Move camera forward */
		break;
	      case SDLK_DOWN:
		/* Move camera backwards */
		break;
	      case SDLK_LEFT:
		/* Move camera left */
		break;
	      case SDLK_RIGHT:
		/* Move camera right */
		break;
	      case SDLK_ESCAPE:
		/* Move camera quit */
		return false;
	      }
	  }  
    }
  return true;
}
