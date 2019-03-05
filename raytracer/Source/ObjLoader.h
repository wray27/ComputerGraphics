#ifndef OBJ_LOADER_H
#define OBJ_LOADER_H


#include <glm/glm.hpp>
#include <vector>
#include "TestModelH.h"
#include <unistd.h>

using namespace std;
using glm::vec2;
using glm::vec3;
using glm::vec4;





// v is a vertex
// vt is the texture coordinate of one vertex
// vn is the normal of one vertex
// f is a face

Triangle scale(Triangle tri){


	float L = 555;
	tri.v0 *= 2/L;
	tri.v1 *= 2/L;
	tri.v2 *= 2/L;

	tri.v0 -= vec4(1,1,1,1);
	tri.v1 -= vec4(1,1,1,1);
	tri.v2 -= vec4(1,1,1,1);

	tri.v0.x *= -1;
	tri.v1.x *= -1;
	tri.v2.x *= -1;

	tri.v0.y *= -1;
	tri.v1.y *= -1;
	tri.v2.y *= -1;

	tri.v0.w = 1.0;
	tri.v1.w = 1.0;
	tri.v2.w = 1.0;
		
	tri.ComputeNormal();
	return tri;
}


bool loadOBJ(const char * path, std::vector<Triangle>&triangles){
	char* buf = (char*)malloc(sizeof(char)*128);
	size_t size = sizeof(char)*128;
	FILE * file = fopen(path, "r");
	char *name = getcwd(buf,size);
	if( file == NULL ){
	    printf("Impossible to open the file !\n");
	    printf("%s\n", name);

	    return false;
	}
	
	

	std::vector< unsigned int > vertexIndices, uvIndices, normalIndices;
	
	std::vector< glm::vec4 > temp_vertices;
	std::vector< glm::vec2 > temp_uvs;
	std::vector< glm::vec4 > temp_normals;
	while( 1 ){

	    char lineHeader[128];
	    // read the first word of the line
	    int res = fscanf(file, "%s", lineHeader);
	    if (res == EOF)
	        break; // EOF = End Of File. Quit the loop.
	    
	    if ( strcmp( lineHeader, "v" ) == 0 ){
	    	// If the first word of the line is “v”, then the rest has to be 3 floats, so create a glm::vec3 out of them, and add it to the vector.
		    glm::vec4 vertex;
		    fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z );
		    vertex.w = 1;
		    temp_vertices.push_back(vertex);
		}
		else if ( strcmp( lineHeader, "vt" ) == 0 ){
			// vt”, then the rest has to be 2 floats, so create a glm::vec2 and add it to the vector.
		    glm::vec2 uv;
		    fscanf(file, "%f %f\n", &uv.x, &uv.y );
		    temp_uvs.push_back(uv);
		}
		else if ( strcmp( lineHeader, "vn" ) == 0 ){
			// for the normal of a vertex
		    glm::vec4 normal;
		    fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z );
		    normal.w = 1;
		    temp_normals.push_back(normal);
		}
		else if ( strcmp( lineHeader, "f" ) == 0 ){
			// for the face
		    std::string vertex1, vertex2, vertex3;
		    unsigned int vertexIndex[3];
		    // unsigned int normalIndex[3];
		    // int matches = fscanf(file, "%d//%d %d//%d %d//%d\n", &vertexIndex[0], &normalIndex[0], &vertexIndex[1], &normalIndex[1], &vertexIndex[2], &normalIndex[2] );
		    int matches = fscanf(file, "%d %d %d\n", &vertexIndex[0], &vertexIndex[1], &vertexIndex[2]);
		    if (matches != 3){
		        printf("File can't be read by our simple parser : ( Try exporting with other options\n");
		        return false;
		    }
		    vertexIndices.push_back(vertexIndex[0]);
		    vertexIndices.push_back(vertexIndex[1]);
		    vertexIndices.push_back(vertexIndex[2]);
		    
		    // uvIndices    .push_back(uvIndex[0]);
		    // uvIndices    .push_back(uvIndex[1]);
		    // uvIndices    .push_back(uvIndex[2]);
		    // normalIndices.push_back(normalIndex[0]);
		    // normalIndices.push_back(normalIndex[1]);
		    // normalIndices.push_back(normalIndex[2]);
		    


		    vec3 white(1.0f,1.0f,1.0f);
			vec4 A;
			vec4 B;
			vec4 C;
			vec4 normal;
			// max and minimum vertex function to find the correct scale!
			for( int i=0; i < vertexIndices.size(); i++ ){
				// For each vertex of each triangle
				unsigned int vertexIndex = vertexIndices[i];
				if(i == 0 ) A = temp_vertices[ vertexIndex-1 ];
				if(i == 1 ) B = temp_vertices[ vertexIndex-1 ];
				if(i == 2 ) C = temp_vertices[ vertexIndex-1 ];
			}


			
			// normal = temp_normals[normalIndex[0]-1];

			Triangle temp_triangle = Triangle( A, B, C, white );
			//temp_triangle = scale(temp_triangle);
			triangles.push_back(temp_triangle);			

		}


	}
	
	return true;

}



#endif
