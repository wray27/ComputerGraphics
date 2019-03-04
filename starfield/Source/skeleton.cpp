#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false


/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
int t;
vector<vec3> stars(1000);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void ColourScreen(screen* screen);
void Draw(screen* screen);
void Interpolate(float a, float b, vector<float>& result);
void Interpolate2(vec3 a, vec3 b, vector<vec3>& result);

int main( int argc, char* argv[] )
{

  
  // float randomNumber;

  for(int i =0; i < stars.size();i++){
    stars[i].x = float (rand())/float(RAND_MAX/2) -1;
    stars[i].y = float (rand())/float(RAND_MAX/2) -1;
    stars[i].z = float (rand())/float(RAND_MAX);


  }
  // for (int i = 0; i < stars.size(); ++i)
  // {
  //   cout << "("
  //   << stars[i].x << "," 
  //   << stars[i].y << "," 
  //   << stars[i].z << ")\n"; //prints the result
  // }
  
  
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  t = SDL_GetTicks();	/*Set start value for timer.*/
  
  while( NoQuitMessageSDL() )
    {
      Draw(screen);
      Update();
      // break;
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}
void Interpolate(float a, float b, vector<float>& result)
{//Most basic version
  for (int i = a; i <= b; ++i)
  {
    result[i-a] = i;
    
  }
}
void Interpolate2(vec3 a, vec3 b, vector<vec3>& result)
{//interpolation with colours 
  float xIncrement = ( b.x - a.x )/(result.size()-1);
  float yIncrement = ( b.y - a.y )/(result.size()-1);
  float zIncrement = ( b.z - a.z )/(result.size()-1);
  for (int i = 0; i <result.size(); ++i)
  {
    result[i].x = a.x + i*xIncrement;
    result[i].y = a.y + i*yIncrement;
    result[i].z = a.z + i*zIncrement;

  }
}
void ColourScreen(screen* screen){
  vec3 red(1.0,0.0,0.0); //top left
  vec3 blue(0.0,0.0,1.0); //top right
  vec3 green(0.0,1.0,0.0); //bottom right
  vec3 yellow(1.0,1.0,0.0); //bottom left
  vec3 magenta(1.0,0.1,0.8);
  vec3 cyan(0.0,0.8,1.0);
  // vec3 colours[6] = {red,blue,green,yellow,magenta,cyan};

  vector<vec3> leftSide(SCREEN_HEIGHT);
  vector<vec3> rightSide(SCREEN_HEIGHT);

  Interpolate2(red,yellow,leftSide);
  Interpolate2(blue,green,rightSide);
  
  // uint32_t x;
  // uint32_t y;

  for(int i=0; i<SCREEN_HEIGHT; i++)
  {
    vector<vec3> row(SCREEN_WIDTH);
    Interpolate2(leftSide[i],rightSide[i],row);
    for (int j = 0; j < SCREEN_WIDTH; ++j)
    { 

      PutPixelSDL(screen, j, i, row[j]);    
    }    
  }

}

/*Place your drawing here*/
void Draw(screen* screen)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  float u[1000];
  float v[1000];
  

  // focal lenght of the camera
  float f = SCREEN_HEIGHT;
  for(size_t s =0; s < stars.size();s++){
    vec3 star= 0.2f * vec3(1,1,1)/(stars[s].z*stars[s].z);
    u[s]= f * (stars[s].x/stars[s].z) +(SCREEN_WIDTH/2);
    v[s]= f * (stars[s].y/stars[s].z) +(SCREEN_HEIGHT/2);
    PutPixelSDL(screen, u[s], v[s], star);

  }
  

  



}

/*Place updates of parameters here*/
void Update()
{
  /* Compute frame time */
  static int t = SDL_GetTicks();
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  float velocity =0.00002;
  for(size_t s =0; s < stars.size();s++){
    stars[s].z = stars[s].z - velocity*dt;
    if(stars[s].z<= 0){
      stars[s].z +=1;
    }
    if(stars[s].z >1){
      stars[s].z -=1;
    }

  }

  
  cout << ".";
  /*Good idea to remove this*/
  // std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/




  // vector<vec3> result(4);
  // vec3 a(1,4,9.2);
  // vec3 b(4,1,9.8);

  // Interpolate2(a,b,result);
  // for (int i = 0; i < result.size(); ++i)
  // {
  //   cout << "("
  //   << result[i].x << "," 
  //   << result[i].y << "," 
  //   << result[i].z << ")";
  // }

  // vector<float> result(10);
  // Interpolate(5,14,result);
  // for (int i = 0; i < result.size(); ++i)
  // {
  //   cout << result[i] << " "; //prints the result
  // }




}

