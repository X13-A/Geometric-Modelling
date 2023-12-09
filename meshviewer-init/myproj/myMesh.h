#pragma once
#include "myFace.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <vector>
#include <string>
#include "myMesh.h"

#include <map>

enum CollapseMode { COLLAPSE_START, COLLAPSE_END, COLLAPSE_AVERAGE };

class myMesh
{
public:
	std::vector<myVertex *> vertices;
	std::vector<myHalfedge *> halfedges;
	std::vector<myFace *> faces;
	std::string name;

	void checkMesh();
	void createHalfEdge(std::map<std::pair<int, int>, myHalfedge*>& twin_map, int prev_i, int i);
	bool readFile(std::string filename);
	void computeNormals();
	void normalize();

	void generateFromCurve(std::vector<myPoint3D> curve, int n, int height);
	void linkCurves(std::vector<myVertex*> vertices1, std::vector<myVertex*> vertices2);

	void subdivisionCatmullClark();

	void splitFaceTRIS(myFace *, myPoint3D *);
	void splitEdge(myHalfedge *, myPoint3D *);
	void splitFaceQUADS(myFace *, myPoint3D *);

	void triangulate();
	bool triangulateFace(myFace *);

	void clear();

	void simplify();
	bool unifyEdge(myHalfedge* he, CollapseMode mode);
	bool unifyEdgeTriangle(myHalfedge* he, CollapseMode mode);
	myHalfedge* getShortestValidEdge();

	void freeVertex(myVertex* v);
	void freePoint(myPoint3D* v);
	void freeHalfEdge(myHalfedge* he);
	void freeFace(myFace* f);
	myMesh(void);
	~myMesh(void);
};

