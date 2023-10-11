#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>
#include <iostream>

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(1.0, 1.0, 1.0);
}

myFace::~myFace(void)
{
	if (normal) delete normal;
}

void myFace::computeNormal()
{
	myVertex* v1 = adjacent_halfedge->source;
	myVertex* v2 = adjacent_halfedge->next->source;
	myVertex* v3 = adjacent_halfedge->prev->source;

	myVector3D u = *v2->point - *v1->point;
	myVector3D v = *v3->point - *v1->point;

	myVector3D n = u.crossproduct(v);

	n.normalize();
	if (normal) delete normal;
	normal = new myVector3D(n);
}
