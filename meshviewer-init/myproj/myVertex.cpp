#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"
#include <iostream>


myVertex::myVertex(void)
{
	point = NULL;
	originof = NULL;
	normal = new myVector3D(1.0,1.0,1.0);
}

myVertex::~myVertex(void)
{
	if (normal) delete normal;
}

void myVertex::computeNormal()
{
	myVector3D n = myVector3D(0, 0, 0);

	// Iterate adjacent faces

	int i = 0;
	myHalfedge* current_he = originof;
	do
	{
		if (!current_he) break;
		if (!current_he->twin) break;
		if (!current_he->adjacent_face) break;
		if (current_he->adjacent_face->normal)
		{
			n += *current_he->adjacent_face->normal;
		}
		current_he = current_he->twin->next;
		i++;
	}
	while (current_he != originof);

	n.normalize();
	if (normal) delete normal;
	normal = new myVector3D(n);
}
