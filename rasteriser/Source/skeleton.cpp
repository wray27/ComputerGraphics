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
void Interpolate( vec2 a, vec2 b, vector<ivec2>& result );
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color );
mat4 setRotationMat(mat4 R,float yaw);
void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen );
void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels );
void testComputePolygonRows();
void DrawPolygonRows( const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels,vec3 color );

void DrawPolygon( const vector<vec4>& vertice,screen* screen,vec3 color );

int main( int argc, char* argv[] )
{
  testComputePolygonRows();
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
    vec3 currentColour;
    
    for(int i=0; i<triangles.size(); i++)
    {
        triangleVerts[0] = triangles[i].v0;
        triangleVerts[1] = triangles[i].v1;
        triangleVerts[2] = triangles[i].v2;
        currentColour = triangles[i].color;

        // DrawPolygonEdges(triangleVerts,screen);
       DrawPolygon(triangleVerts,screen,currentColour);
        
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
void ComputePolygonRows(
       const vector<ivec2>& vertexPixels,
       vector<ivec2>& leftPixels,
       vector<ivec2>& rightPixels )
{


    
    // 1. Find max and min y-value of the polygon
    //    and compute the number of rows it occupies.
    int largest =  -numeric_limits<int>::max();
    int smallest= +numeric_limits<int>::max();
    for(int i = 0; i < vertexPixels.size(); i ++){
        if(vertexPixels[i].y < smallest) smallest = vertexPixels[i].y;
        if(vertexPixels[i].y > largest) largest = vertexPixels[i].y;
    }
    
    int ROWS = largest - smallest +1;
    // cout << ROWS <<endl;



    // 2. Resize leftPixels and rightPixels
    //    so that they have an element for each row.
    leftPixels.resize(ROWS);
    rightPixels.resize(ROWS);
    // cout << leftPixels.size() <<endl;
    



  
    // 3. Initialize the x-coordinates in leftPixels
    //    to some really large value and the x-coordinates
    //    in rightPixels to some really small value.
    
    // initializes left and right pixels row to contain the smallest and largest
    // x values of the pixels in that row respectively
    
    for( int i=0; i<ROWS; ++i )
    {
       leftPixels[i].x  = +numeric_limits<int>::max();
       rightPixels[i].x = -numeric_limits<int>::max();
    }

    // 4. Loop through all edges of the polygon and use
    //    linear interpolation to find the x-coordinate for
    //    each row it occupies. Update the corresponding
    // values in rightPixels and leftPixels. 

    ivec2 delta;
    
    ivec2 a;
    ivec2 b;
    int V = vertexPixels.size();
    // vector< vector<ivec2> > edges(V);
    
    
    for( int i=0; i<V; ++i )
    {
        int j = (i+1)%V; // The next vertex
        a = vertexPixels[i];
        b = vertexPixels[j];

        delta.x = glm::abs( a.x - b.x );
        delta.y = glm::abs( a.y - b.y );
        int pixels = glm::max( delta.x, delta.y ) + 1;
        
        // edges[i].resize(pixels);
        
        vector<ivec2> line( pixels );
        Interpolate( a, b, line );


        for(int p = 0; p < pixels;p++){
            for(int j =0 ; j < ROWS ; j++){
                if(line[p].x < leftPixels[j].x && line[p].y == j + smallest){
                   leftPixels[j].x = line[p].x;
                   

                   leftPixels[j].y = smallest + j;
                } 
                

                if(line[p].x > rightPixels[j].x && line[p].y == j + smallest){
                    rightPixels[j].x = line[p].x;
                    

                    rightPixels[j].y = smallest + j;
                } 


            }

        }
        // edges[i] = line;    
    } 

    
    // for(int i = 0; i < edges.size();i ++){
    //     for(int j =0 ; j < ROWS; j++){
    //         if(edges[i][j].x < leftPixels[j].x && edges[i][j].y == j + smallest){
    //            leftPixels[j].x = edges[i][j].x;
               

    //            leftPixels[j].y = smallest + j;
    //         } 
            

    //         if(edges[i][j].x > rightPixels[j].x && edges[i][j].y == j + smallest){
    //             rightPixels[j].x = edges[i][j].x;
                

    //             rightPixels[j].y = smallest + j;
    //         } 
    //     }
    // }
    // cout << edges.size() <<endl;


}
void testComputePolygonRows(){
    vector<ivec2> vertexPixels(3);
    vertexPixels[0] = ivec2(10, 5);
    vertexPixels[1] = ivec2( 5,10);
    vertexPixels[2] = ivec2(15,15);
    vector<ivec2> leftPixels;
    vector<ivec2> rightPixels;
    ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
    // cout << leftPixels.size() <<endl;
    for( int row=0; row<leftPixels.size(); ++row )
    {
        
        cout << "Start: ("
            << leftPixels[row].x << ","
            << leftPixels[row].y << "). "
            << "End: ("
            << rightPixels[row].x << ","
            << rightPixels[row].y << "). " << endl;
    }
}

void Interpolate( vec2 a, vec2 b, vector<ivec2>& result )
{
    int N = result.size();
    float temp = float(max(N-1,1));
    
    vec2 step;

    step.x = (b.x - a.x) / temp;
    step.y = (b.y - a.y) / temp;

    vec2 current( a );
    for( int i=0; i<N; ++i )
    {
       result[i].x = current.x;
       result[i].y = current.y;
       current.x += step.x;

       current.y += step.y;
       // cout << result[i].x << endl;
       // cout << result[i].y << endl;
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
void DrawPolygonRows( const vector<ivec2>& leftPixels,
                  const vector<ivec2>& rightPixels, vec3 color,screen* screen ){
    for(int i = 0; i < leftPixels.size(); i ++){
        int start = leftPixels[i].x;
        int stop =  rightPixels[i].x;
        for(int j = start; j <= stop; j++){
            PutPixelSDL( screen, start + j , i, color );
        }
    }

}

void DrawPolygon( const vector<vec4>& vertices,screen* screen,vec3 color)
{
       int V = vertices.size();
       vector<ivec2> vertexPixels( V );
       for( int i=0; i<V; ++i )
           VertexShader( vertices[i], vertexPixels[i] );
       vector<ivec2> leftPixels;
       vector<ivec2> rightPixels;
       ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
       DrawPolygonRows( leftPixels, rightPixels,color,screen);
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