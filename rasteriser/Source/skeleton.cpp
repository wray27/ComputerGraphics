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
#define PI 3.14159265358979323546

vector<Triangle> triangles;
vec4 cameraPos( 0, 0, -3.001,1 );
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
float yaw = 0;
vec4 lightPos(0,-0.5,-0.7,1);
vec3 lightPower = 18.0f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );
struct Pixel
{
       int x;
       int y;
       float zinv;
       vec4 pos3d;
};


struct Vertex
{
    vec4 position;
    bool onScreen;

    // vec4 normal;
    // vec3 reflectance;
};
struct ClippedTriangle
{
       vector<Vertex> vertices;
       int offScreenCount;
};

vec4 currentNormal;
vec3 currentReflectance;
// focal length of the camera
float f = SCREEN_WIDTH;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& v, ivec2& p );
void VertexShader( const Vertex& v, Pixel& p);
void Interpolate( vec2 a, vec2 b, vector<ivec2>& result );
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color );
mat4 setRotationMat(mat4 R,float yaw);
void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen );
void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels );
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels );
void testComputePolygonRows();
void DrawPolygonRows( const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels,vec3 color );
void DrawPolygonRows( const vector<Pixel>& leftPixels,const vector<Pixel>& rightPixels,screen* screen );
void PixelShader( const Pixel& p,screen* screen);
vec3 light( const vec4& v );
void DrawPolygon( const vector<Vertex>& vertices,screen* screen);
vec4 toClipSpace(vec4 worldSpace);
bool onScreen(vec4 clipSpace);

int main( int argc, char* argv[] ){
  // testComputePolygonRows();
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
bool onScreen(vec4 clipSpace){
    bool onScreen = false;
    
    if( clipSpace.x >= clipSpace.w*(SCREEN_WIDTH/2)*-1 &&  
        clipSpace.x <=  clipSpace.w*(SCREEN_WIDTH/2)  &&
        clipSpace.y >= clipSpace.w*(SCREEN_HEIGHT/2)*-1  &&  
        clipSpace.y <= clipSpace.w*(SCREEN_HEIGHT/2) ){
        onScreen = true;
    }
    return onScreen;
}
// Changes the coordinate from world space to clip space
vec4 toClipSpace(vec4 worldSpace){
    vec4 clipSpace;
    clipSpace = worldSpace - cameraPos;
    clipSpace.w = clipSpace.z/f;
    // clipSpace -= cameraPos;
    return clipSpace;

}
/*Place your drawing here*/
void Draw(screen* screen){
    /* Clear buffer */
    memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
    
    
    // clears the depth buffer
    for( int y=0; y<SCREEN_HEIGHT; ++y )
        for( int x=0; x<SCREEN_WIDTH; ++x )
            depthBuffer[y][x] = 0;
   
    vec3 colour(1.0,1.0,1.0);
    vector<Vertex> triangleVerts(3);
    vector<vec4> clipPositions(3);
    
    vector<Triangle> keepTriangles;
    vector<ClippedTriangle> clippedTriangles;

    

    
    // ivec2 delta = glm::abs( a - b );
    // int pixels = glm::max( delta.x, delta.y ) + 1;
    
    // vector<ivec2> line( pixels );
    // Interpolate( a, b, line );

    // boolean to add triangles that do not need to be clipped
    bool addTriangle;
    int offScreenCount = 0;
    
    // CLIPPING
    for(int i=0; i<triangles.size(); i++){
       

        addTriangle = true;
        triangleVerts[0].position = triangles[i].v0;
        triangleVerts[1].position = triangles[i].v1;
        triangleVerts[2].position = triangles[i].v2;
        //converts the coordinates of the vertices of each triangle from world space to clip space
        for(int j = 0; j <3;j++){
            clipPositions[j] = toClipSpace(triangleVerts[j].position); 
            if(!onScreen(clipPositions[j])){
                triangleVerts[j].onScreen = false;
                addTriangle = false;
                offScreenCount++;
            } else{
                triangleVerts[j].onScreen = true;
            }
        }
        
        
        if(addTriangle){
            keepTriangles.push_back(triangles[i]);
        }else{
            ClippedTriangle clippedTri;
            clippedTri.vertices = triangleVerts;
            clippedTri.offScreenCount = offScreenCount;

            clippedTriangles.push_back(clippedTri);
        }

        offScreenCount = 0;

    }

    vector<vec4> intersections(2);
    
    // ADDING NEW TRIANGLES
    for(int i = 0; i < clippedTriangles.size();i++){
        
    // CASE 1 
    /*
        If there's 1 vertex out, then you gotta take the 2 intersection points and 2 original vertices and create 2 new triangles

    */
        if(clippedTriangles[i].offScreenCount == 1){
            int offScreenIndex = 0;
            for(int j = 0; j < 3;j++) if(!clippedTriangles[i].vertices[j].onScreen) offScreenIndex = j;;
            for(int j = 0; j < 3;j++){
                if(j != offScreenIndex){
                    // clippedTriangles[i].vertices[j]
                }
            }
        }



    // CASE 2 
    /*
        If 2 vertices are out, just join the intersection point to create new triangle. 

    */

        if(clippedTriangles[i].offScreenCount == 2){
            int onScreenIndex = 0;
            for(int j = 0; j < 3;j++) if(clippedTriangles[i].vertices[j].onScreen) onScreenIndex = j;;
            for(int j = 0; j < 3;j++){
                if(j == onScreenIndex){
                    // clippedTriangles[i].vertices[j]
                }
            }
        }

        
    }
    





    for(int i=0; i<keepTriangles.size(); i++)
    {  
        currentReflectance = keepTriangles[i].color;
        currentNormal = keepTriangles[i].normal;
        
        triangleVerts[0].position = keepTriangles[i].v0;
        triangleVerts[1].position = keepTriangles[i].v1;
        triangleVerts[2].position = keepTriangles[i].v2;

        // power = light(triangleVerts[0]);
        // reflectance = (power + indirectLightPowerPerArea) * currentColor;
        // triangleVerts[0].reflectance = currentColor ;
        // triangleVerts[0].normal = triangles[i].normal;

        
        // power = light(triangleVerts[1]);
        // reflectance = (power + indirectLightPowerPerArea) * currentColor;
        // triangleVerts[1].reflectance = currentColor ;
        // triangleVerts[1].normal = triangles[i].normal;


        // power = light(triangleVerts[2]);
        // reflectance = (power + indirectLightPowerPerArea) * currentColor;
        // triangleVerts[2].reflectance = currentColor;
        // triangleVerts[2].normal = triangles[i].normal;

        

        // DrawPolygonEdges(triangleVerts,screen);
        DrawPolygon(triangleVerts,screen);
        
    }
}
void PixelShader( const Pixel& p,screen* screen){
    int x = p.x;
    int y = p.y;

    vec3 power = light(p.pos3d);
    vec3 illumination = (power + indirectLightPowerPerArea) * currentReflectance;
    // cout << "colour val  = " << p.illumination.x;
    if( p.zinv > depthBuffer[y][x] )
    {

        depthBuffer[y][x] = p.zinv;
        PutPixelSDL( screen, x, y, illumination);
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
vec3 light( const vec4& v ){
    vec3 colour = vec3(0.0f,0.0f,0.0f);

    //calculatiing the distance from light source - r
    float distance = glm::distance(lightPos,v);
    
    // unit vector describing the direction from the surface point to the light source
    vec4 _r = glm::normalize(lightPos - v);
    // vec3 r = glm::normalize(vec3(_r.x,_r.y,_r.z));

    //normal pointing out of the surface as a unit vector
    vec4 _n = glm::normalize(currentNormal);
    // vec3 n = glm::normalize(vec3(triangles[i.triangleIndex].normal.x,triangles[i.triangleIndex].normal.y,triangles[i.triangleIndex].normal.x));


    // Adds light to the scene withno colour
    colour = (vec3)( lightPower * max(glm::dot(_r,_n),0.0f) ) / (float)(4.0f*PI* pow((double)distance,2));
    // if ( distance <=  0 ) cout <<  "what is going on";
    // cout << " colour : " << colour.y;
    // adds colour tot the scene
    // colour = colour * triangles[i.triangleIndex].color;


    return colour;
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
void VertexShader( const Vertex& v, Pixel& p){
    vec4 temp;
    glm::mat4 R;

    p.pos3d.x = v.position.x;
    p.pos3d.y = v.position.y;
    p.pos3d.z = v.position.z;
    p.pos3d.w = v.position.w;

    // vec3 power = light(v);
    // p.illumination = (power + indirectLightPowerPerArea) * v.reflectance;
    // // cout << "colour val  = " << p.illumination.x;
    // p.illumination = vec3(1.0f,1.0f,1.0f);
   
    
    


    R = setRotationMat(R,yaw);
    temp = R*v.position;

    // puts the coordinates of the vertex in the same system as the camera
    temp.x = temp.x - cameraPos.x;
    temp.y = temp.y - cameraPos.y;
    temp.z = temp.z - cameraPos.z;
    // calculates the inverse of the z coordinate for each pixel - used for depth calculations
    p.zinv = 1 / temp.z;

    


    p.x = (f * (temp.x/temp.z)) + (SCREEN_WIDTH /2);
    p.y = (f * (temp.y/temp.z)) + (SCREEN_HEIGHT /2);
 
}
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels ){

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
    
    Pixel a;
    Pixel b;
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
        
        vector<Pixel> line( pixels );
        Interpolate( a, b, line );


        for(int p = 0; p < pixels;p++){
            for(int j =0 ; j < ROWS ; j++){
                
                if(line[p].x < leftPixels[j].x && line[p].y == j + smallest){
                   leftPixels[j].x = line[p].x;
                   leftPixels[j].zinv = line[p].zinv;
                   leftPixels[j].y = smallest + j;
                   leftPixels[j].pos3d = line[p].pos3d;
                } 

                if(line[p].x > rightPixels[j].x && line[p].y == j + smallest){
                    rightPixels[j].x = line[p].x;
                    rightPixels[j].zinv = line[p].zinv;
                    rightPixels[j].y = smallest + j;
                    rightPixels[j].pos3d = line[p].pos3d;
                } 


            }

        }
        // edges[i] = line;    
    } 
}
void ComputePolygonRows(const vector<ivec2>& vertexPixels,vector<ivec2>& leftPixels,  vector<ivec2>& rightPixels ){

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
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result ){
    Pixel c;
    Pixel d;
    int N = result.size();
    float temp = float(max(N-1,1));
    
    vec3 step;
    vec4 posStep;
    c.x = a.x;
    c.y = a.y;
    c.zinv = a.zinv;
    c.pos3d = a.pos3d;

    d.x = b.x;
    d.y = b.y;
    d.zinv = b.zinv;
    d.pos3d = b.pos3d;



    c.pos3d.x = c.pos3d.x * c.zinv;
    c.pos3d.y  = c.pos3d.y * c.zinv;
    
    d.pos3d.x = d.pos3d.x * d.zinv;
    d.pos3d.y = d.pos3d.y * d.zinv;



    step.x = (b.x - a.x) / temp;
    step.y = (b.y - a.y) / temp;
    step.z = (b.zinv - a.zinv) / temp;

    posStep.x = (d.pos3d.x - c.pos3d.x) / temp;
    posStep.y = (d.pos3d.y - c.pos3d.y) / temp;
    posStep.z = (b.pos3d.z - a.pos3d.z) / temp;



    vec3 current( a.x,a.y,a.zinv );
    vec4 currentPos(c.pos3d.x,c.pos3d.y,c.pos3d.z,1.0f);
    for( int i=0; i<N; ++i )
    {
       result[i].x = current.x;
       result[i].y = current.y;
       result[i].zinv = current.z;

       result[i].pos3d.x = currentPos.x;
       result[i].pos3d.y = currentPos.y;
       result[i].pos3d.z = currentPos.z;


       
       current.x += step.x;
       current.y += step.y;
       current.z += step.z;

       currentPos.x += posStep.x;
       currentPos.y += posStep.y;
       currentPos.z += posStep.z;
       currentPos.w = 0;
       // cout << result[i].x << endl;
       // cout << result[i].y << endl;
    }
    for(int i = 0; i < N;i ++)
    {
        result[i].pos3d.x = result[i].pos3d.x / result[i].zinv;
        result[i].pos3d.y = result[i].pos3d.y / result[i].zinv;
    }

}
void Interpolate( vec2 a, vec2 b, vector<ivec2>& result ){
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
void DrawPolygonRows( const vector<ivec2>& leftPixels,const vector<ivec2>& rightPixels, vec3 color,screen* screen ){
    for(int i = 0; i < leftPixels.size(); i ++){
        int start = leftPixels[i].x;
        int stop =  rightPixels[i].x;


        for(int j = start; j <= stop; j++){
            PutPixelSDL( screen, j , leftPixels[i].y, color );
        }
    }
}
void DrawPolygonRows( const vector<Pixel>& leftPixels,const vector<Pixel>& rightPixels,screen* screen ){
    int start;
    int stop;
  
    for(int i = 0; i < leftPixels.size(); i ++){
        start = leftPixels[i].x;
        stop = rightPixels[i].x;
        vector<Pixel> currrentRow(stop - start + 1);
        Interpolate(leftPixels[i],rightPixels[i], currrentRow);
        for(int j = 0; j < currrentRow.size(); j++){
            // exits early if pixels are no longer on the screen
            if( currrentRow[j].x < 0 || currrentRow[j].x >= SCREEN_WIDTH||leftPixels[i].y < 0 || leftPixels[i].y >=SCREEN_HEIGHT) continue;
            // A new pixel is drawn if the inverse of the z coordinate is larger than the value in the depth buffer
            // if( currrentRow[j].zinv >= depthBuffer[currrentRow[j].x][leftPixels[i].y] ){
                
            //     PutPixelSDL( screen, currrentRow[j].x,currrentRow[j].y, color );
            //     // stores new pixel that is drawn in depth buffer
            //     depthBuffer[currrentRow[j].x][currrentRow[j].y] = currrentRow[j].zinv;

            // }  
            PixelShader(currrentRow[j],screen);      
        }
    }
}
void DrawPolygon( const vector<Vertex>& vertices,screen* screen){
       int V = vertices.size();
       vector<Pixel> vertexPixels( V );
       for( int i=0; i<V; ++i )
           VertexShader( vertices[i], vertexPixels[i]);
      
       vector<Pixel> leftPixels;
       vector<Pixel> rightPixels;
       ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
       DrawPolygonRows( leftPixels, rightPixels,screen);
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
bool Update(){
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
                    case SDLK_w:
                    {
                        lightPos += forward;
                        break;
                    }
                    case SDLK_s:
                    {
                        lightPos -= forward;
                        break;
                    }
                    case SDLK_a:
                    {
                        lightPos += right;
                        break;
                    }
                    case SDLK_d:
                    {
                        lightPos -= right;
                        break;
                    }
                    case SDLK_e:
                    {
                        lightPos -= down;
                        break;
                    }
                    case SDLK_q:
                    {
                        lightPos += down;
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