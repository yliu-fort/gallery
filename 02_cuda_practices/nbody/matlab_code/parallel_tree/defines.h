#ifndef DEFINE_H
#define DEFINE_H

#include <stdint.h>

typedef struct
{
   double x,y,z,w;
} double4;

typedef struct
{
   float x,y,z,w;
} float4;

typedef struct
{
	float4 bbox; // store min coordinates and scale for bounding box
	int begin, end;
	uint32_t level;

	//Node* childA;
	//Node* childB;
    int childA;
    int childB;
    
	//uint32_t mask; 
    bool isOctreeNode;
	bool isLeafNode;
    
    float4 scalar;
} Node;



#endif