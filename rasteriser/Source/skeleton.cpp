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
#define PI 3.14159255358979323546

vector<Triangle> triangles;
vec4 cameraPos( 0, 0, -3.001,1 );
float yaw = 0;
// focal length of the camera
float f = SCREEN_WIDTH;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& v, ivec2& p );
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color );
mat4 setRotationMat(mat4 R,float yaw);
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
    glm::mat4 R;
    R = setRotationMat(R,yaw);
    temp = R*v;
    temp.x = temp.x - cameraPos.x;
    temp.y = temp.y - cameraPos.y;
    temp.z = temp.z - cameraPos.z;



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
    static int t = SDL_GetTicks();
    /* Compute frame time */
    int t2 = SDL_GetTicks();
    // float dt = float(t2-t);
    t = t2;
    float change =0.1f;

    mat4 R;
    R = setRotationMat(R,yaw);
    vec4 forward( R[2][0], R[2][1], R[2][2], 1 );
    vec4 right(   R[0][0], R[0][1], R[0][2], 1 );
    vec4 down(    R[1][0], R[1][1], R[1][2], 1 );


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
                    // case SDLK_w:
                    // {
                    //     lightPos += forward;
                    //     break;
                    // }
                    // case SDLK_s:
                    // {
                    //     lightPos -= forward;
                    //     break;
                    // }
                    // case SDLK_a:
                    // {
                    //     lightPos += right;
                    //     break;
                    // }
                    // case SDLK_d:
                    // {
                    //     lightPos -= right;
                    //     break;
                    // }
                    // case SDLK_e:
                    // {
                    //     lightPos -= down;
                    //     break;
                    // }
                    // case SDLK_q:
                    // {
                    //     lightPos += down;
                    //     break;
                    // }
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