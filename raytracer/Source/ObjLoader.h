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



// not 100% complete, need to ensure vt and vn's are pushed onto the Triangles vector


// v is a vertex
// vt is the texture coordinate of one vertex
// vn is the normal of one vertex
// f is a face

std::vector<glm::vec4> scale(std::vector<glm::vec4>verticies) {

	float max_x = -10.0;
	float max_y = -10.0;
	float max_z = -10.0;

	float min_x = 10.0;
	float min_y = 10.0;
	float min_z = 10.0;

	std::vector<glm::vec4> output;	

	for(int i = 0; i < verticies.size(); i++) {

		if(verticies[i].x > max_x) {
			max_x = verticies[i].x;
		}
		if(verticies[i].y > max_y) {
			max_y = verticies[i].y;
		}
		if(verticies[i].z > max_z) {
			max_z = verticies[i].z;
		}

		if(verticies[i].x < min_x) {
			min_x = verticies[i].x;
		}

		if(verticies[i].x < min_y) {
			min_y = verticies[i].y;
		}

		if(verticies[i].x < min_z) {
			min_z = verticies[i].z;
		}
	}

	float max_diff_x = max_x-min_x;
	float max_diff_y = max_y-min_y;
	float max_diff_z = max_z-min_z;

	float max = std::max(max_diff_x, max_diff_y);
	max = std::max(max, max_diff_z);

	float scale = (0.7 / max) - 1;

	for(int i = 0; i < verticies.size(); i++) {
		output.push_back((verticies[i] * scale) + vec4(0.35,0.5,-0.5,0));

	}

	cout << "Here" << endl;

	return output;
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

	std::vector<Triangle> scaledTriangles;
	

	std::vector< unsigned int > vertexIndicesA, vertexIndicesB,vertexIndicesC, uvIndices, normalIndices;
	
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
		    vertex.y = -vertex.y;
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
		    vertexIndicesA.push_back(vertexIndex[0]);
		    vertexIndicesB.push_back(vertexIndex[1]);
		    vertexIndicesC.push_back(vertexIndex[2]);
		    
		    // uvIndices    .push_back(uvIndex[0]);
		    // uvIndices    .push_back(uvIndex[1]);
		    // uvIndices    .push_back(uvIndex[2]);
		    // normalIndices.push_back(normalIndex[0]);
		    // normalIndices.push_back(normalIndex[1]);
		    // normalIndices.push_back(normalIndex[2]);
		  

			// for( int i=0; i < vertexIndices.size(); i++ ){
			// 	// For each vertex of each triangle
			// 	// cout << vertexIndices[i] << endl;
			// 	unsigned int vertexIndex = vertexIndices[i];
			// 	if(i == 0 ) A = temp_vertices[ vertexIndex-1 ];
			// 	if(i == 1 ) B = temp_vertices[ vertexIndex-1 ];
			// 	if(i == 2 ) C = temp_vertices[ vertexIndex-1 ];
			// }


			// // normal = temp_normals[normalIndex[0]-1];



			// Triangle temp_triangle = Triangle( A, B, C, white );
			// //temp_triangle = scale(temp_triangle);
			// //triangles.push_back(temp_triangle);
			// scaledTriangles.push_back(temp_triangle);

		}

		


	}

	vec3 white(0.5f,0.5f,0.5f);
	vec4 A;
	vec4 B;
	vec4 C;
	vec4 normal;

	std::vector<glm::vec4> scaled_verticies;
	//scale temp_verticies:
	scaled_verticies = scale(temp_vertices);
	cout << scaled_verticies[0].x << endl;

	for(int i = 0; i < vertexIndicesA.size(); i++){
		A = scaled_verticies[vertexIndicesA[i]-1];
		B = scaled_verticies[vertexIndicesB[i]-1];
		C = scaled_verticies[vertexIndicesC[i]-1];

		Triangle temp_triangle = Triangle(A, B, C, white,false);
		triangles.push_back(temp_triangle);

	}


	// cout << scaledTriangles[1].v0.x<< "," << scaledTriangles[1].v0.y << "," << scaledTriangles[1].v0.z <<endl;	
	return true;

}



#endif
