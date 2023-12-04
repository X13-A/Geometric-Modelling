#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myvector3d.h"
#include <stdlib.h>
#include <set>
#include <iostream>
#include <algorithm>

using namespace std;

myMesh::myMesh(void)
{
	/**** TODO ****/
}


myMesh::~myMesh(void)
{
	/**** TODO ****/
}

void myMesh::clear()
{
	for (unsigned int i = 0; i < vertices.size(); i++) if (vertices[i]) delete vertices[i];
	for (unsigned int i = 0; i < halfedges.size(); i++) if (halfedges[i]) delete halfedges[i];
	for (unsigned int i = 0; i < faces.size(); i++) if (faces[i]) delete faces[i];

	vector<myVertex *> empty_vertices;    vertices.swap(empty_vertices);
	vector<myHalfedge *> empty_halfedges; halfedges.swap(empty_halfedges);
	vector<myFace *> empty_faces;         faces.swap(empty_faces);
}

void myMesh::checkMesh()
{
	int null_next = 0;
	int null_prev = 0;
	int null_twin = 0;
	int null_face = 0;
	int null_source = 0;
	int incoherent_next = 0;
	int incoherent_prev = 0;
	int incoherent_twin = 0;

	for (myHalfedge* he : halfedges)
	{
		if (he->twin == nullptr) null_twin++;
		if (he->next == nullptr) null_next++;
		if (he->prev == nullptr) null_prev++;
		if (he->source == nullptr) null_source++;
		if (he->adjacent_face == nullptr) null_face++;
		if (he->next->prev != he) incoherent_next++;
		if (he->prev->next != he) incoherent_prev++;
		if (he->twin && he->twin->twin != he) incoherent_twin++;
	}

	int null_originof = 0;
	int incoherent_originof = 0;

	for (myVertex* v : vertices)
	{
		if (v->originof == nullptr) null_originof++;
		if (v->originof != nullptr && v->originof->source != v) incoherent_originof++;
	}

	int incoherent_face_next_prev = 0;
	for (myFace* face : faces)
	{
		int nexts = 0;
		myHalfedge* he = face->adjacent_halfedge;
		do
		{
			he = he->next;
			nexts++;
		} while (he != face->adjacent_halfedge);

		int prevs = 0;
		he = face->adjacent_halfedge;
		do
		{
			he = he->prev;
			prevs++;
		} while (he != face->adjacent_halfedge);

		if (nexts != prevs) incoherent_face_next_prev++;
	}

	// Check if sources are ok
	int duplicate_source_in_triangle_count = 0;
	for (myHalfedge* he : halfedges)
	{
		if (he->source == he->next->source)
		{
			duplicate_source_in_triangle_count++;
		}
	}

	// Check if all faces are closed
	int he_count = 0;
	for (myFace* face : faces)
	{
		myHalfedge* he = face->adjacent_halfedge;
		do
		{
			he_count++;
			he = he->next;
		} while (he != face->adjacent_halfedge);
	}

	int he_delta = he_count - halfedges.size();
	if (he_delta != 0)
	{
		std::cout << "Warning! Incoherent half-edge number (" << he_delta << ")" << std::endl;
	}
	if (duplicate_source_in_triangle_count != 0)
	{
		std::cerr << "Warning! Some halfedges have the same source than their next (" << duplicate_source_in_triangle_count << ")!" << std::endl;
	}
	if (null_next != 0)
	{
		std::cout << "Warning! Not all edges have their next (" << null_next << ")!\n";
	}
	if (null_prev != 0)
	{
		std::cout << "Warning! Not all edges have their prev (" << null_prev << ")!\n";
	}
	if (null_twin != 0)
	{
		std::cout << "Warning! Not all edges have their twin (" << null_twin << ")!\n";
	}
	if (null_face != 0)
	{
		std::cout << "Warning! Not all edges have their face (" << null_face << ")!\n";
	}
	if (null_source != 0)
	{
		std::cout << "Warning! Not all edges have their source (" << null_source << ")!\n";
	}
	if (incoherent_next != 0)
	{
		std::cout << "Warning! Incoherent next (" << incoherent_next << ")!\n";
	}
	if (incoherent_prev != 0)
	{
		std::cout << "Warning! Incoherent prev (" << incoherent_prev << ")!\n";
	}
	if (incoherent_twin != 0)
	{
		std::cout << "Warning! Incoherent twin (" << incoherent_twin << ")!\n";
	}
	if (null_originof != 0)
	{
		std::cout << "Warning! Not all vertices have their originof (" << null_originof << ")!\n";
	}
	if (incoherent_originof != 0)
	{
		std::cout << "Warning! Incoherent originof (" << incoherent_originof << ")!\n";
	}
	if (incoherent_face_next_prev != 0)
	{
		std::cout << "Warning! Some faces's nexts and prevs don't make sense (" << incoherent_face_next_prev << ")!\n";
	}
	if (null_next + null_prev + null_twin + null_face + null_source + incoherent_next + incoherent_prev + incoherent_twin + null_originof + incoherent_originof + incoherent_face_next_prev == 0) std::cout << "Each edge is set!\n";
}

float degreesToRadians(float degrees) 
{
	return degrees * (3.14159265359 / 180.0);
}

myPoint3D getRotatedPoint(const myPoint3D& origin, const myPoint3D& p, float theta) 
{
	// Step 1: Translate point p to the origin
	float translatedX = p.X - origin.X;
	float translatedZ = p.Z - origin.Z;

	// Step 2: Rotate the point around the Y axis
	float rotatedX = translatedX * cos(theta) - translatedZ * sin(theta);
	float rotatedZ = translatedX * sin(theta) + translatedZ * cos(theta);

	// Step 3: Translate the point back
	myPoint3D rotatedPoint;
	rotatedPoint.X = rotatedX + origin.X;
	rotatedPoint.Y = p.Y; // y coordinate remains unchanged
	rotatedPoint.Z = rotatedZ + origin.Z;

	return rotatedPoint;
}


std::vector<myVertex*> createVerticesFromPoints(std::vector<myPoint3D> points)
{
	std::vector<myVertex*> vertices;
	for (int i = 0; i < points.size(); i++)
	{
		myVertex* v = new myVertex();
		v->point = new myPoint3D(points.at(i));
		vertices.push_back(v);
	}
	return vertices;
}


double clamp(double value, double min, double max)
{
	return std::max(std::min(value, max), min);
}

double smoothstep(double edge0, double edge1, double x) 
{
	x = clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f);
	return x * x * (3 - 2 * x);
}

void myMesh::generateFromCurve()
{
	int height = 100;
	int n = 50;
	if (n <= 0) return;

	std::vector<myPoint3D> curve;
	for (size_t i = 0; i < height; i++)
	{
		float x = 0.25f + ((sin(float(i)/4) + 1)/2)/8;
		float y = -0.5 + float(i)/ height;
		float z = 0;
		curve.push_back(myPoint3D(x, y ,z));
	}

	// Create first curve vertices
	std::vector<std::vector<myVertex*>> curves_vertices = { createVerticesFromPoints(curve) };

	for (int i = 1; i < n; i++)
	{
		// Get new curve by rotating the original
		std::vector<myPoint3D> next_curve;
		for (int j = 0; j < curve.size(); j++)
		{
			next_curve.push_back(getRotatedPoint(myPoint3D(0, 0, 0), curve[j], degreesToRadians(360.0) * i / n));
		}

		std::vector<myVertex*> current_curve_vertices = curves_vertices.back();
		std::vector<myVertex*> next_curve_vertices = createVerticesFromPoints(next_curve);

		linkCurves(current_curve_vertices, next_curve_vertices);
		curves_vertices.push_back(next_curve_vertices);
	}
	linkCurves(curves_vertices.back(), curves_vertices[0]);
}

void myMesh::linkCurves(std::vector<myVertex*> vertices1, std::vector<myVertex*> vertices2)
{
	// Set shortest and longest vector
	std::vector<myVertex*> shortest_vector;
	std::vector<myVertex*> longest_vector;

	if (vertices1.size() > vertices2.size())
	{
		longest_vector = vertices1;
		shortest_vector = vertices2;
	}
	else
	{
		longest_vector = vertices2;
		shortest_vector = vertices1;
	}

	// Triangulate
	myHalfedge* last_boundary_edge = nullptr;
	myVertex* last_shortest_vector_vertex = shortest_vector[shortest_vector.size() - 1];
	for (int i = 0; i < longest_vector.size(); i++)
	{
		if (i < shortest_vector.size() - 1)
		{
			// Create first triangle
			myHalfedge* he1 = new myHalfedge;
			myHalfedge* he2 = new myHalfedge;
			myHalfedge* he3 = new myHalfedge;

			he1->source = shortest_vector[i];
			shortest_vector[i]->originof = he1;

			he2->source = longest_vector[i + 1];
			longest_vector[i + 1]->originof = he2;

			he3->source = longest_vector[i];
			longest_vector[i]->originof = he3;

			he1->next = he2;
			he2->prev = he1;

			he2->next = he3;
			he3->prev = he2;

			he3->next = he1;
			he1->prev = he3;

			// Set first face
			myFace* face1 = new myFace();
			face1->adjacent_halfedge = he1;
			he1->adjacent_face = face1;
			he2->adjacent_face = face1;
			he3->adjacent_face = face1;

			// Create second triangle
			myHalfedge* he4 = new myHalfedge;
			myHalfedge* he5 = new myHalfedge;
			myHalfedge* he6 = new myHalfedge;

			he4->source = shortest_vector[i];
			shortest_vector[i]->originof = he4;

			he5->source = shortest_vector[i + 1];
			shortest_vector[i + 1]->originof = he5;

			he6->source = longest_vector[i + 1];
			longest_vector[i + 1]->originof = he6;

			he4->next = he5;
			he5->prev = he4;

			he5->next = he6;
			he6->prev = he5;

			he6->next = he4;
			he4->prev = he6;

			// Set second face
			myFace* face2 = new myFace();
			face2->adjacent_halfedge = he4;
			he4->adjacent_face = face2;
			he5->adjacent_face = face2;
			he6->adjacent_face = face2;

			// Set twins
			he1->twin = he6;
			he6->twin = he1;

			if (last_boundary_edge != nullptr)
			{
				he4->twin = last_boundary_edge;
				last_boundary_edge->twin = he4;
			}
			last_boundary_edge = he4;

			// Store new edges and faces
			faces.push_back(face1);
			faces.push_back(face2);
			halfedges.push_back(he1);
			halfedges.push_back(he2);
			halfedges.push_back(he3);
			halfedges.push_back(he4);
			halfedges.push_back(he5);
			halfedges.push_back(he6);
		}
		else if (i < longest_vector.size() - 1)
		{
			// Create triangle
			myHalfedge* he1 = new myHalfedge;
			myHalfedge* he2 = new myHalfedge;
			myHalfedge* he3 = new myHalfedge;

			he1->source = last_shortest_vector_vertex;
			last_shortest_vector_vertex->originof = he1;

			he2->source = longest_vector[i + 1];
			longest_vector[i + 1]->originof = he2;

			he3->source = longest_vector[i];
			longest_vector[i]->originof = he3;

			he1->next = he2;
			he2->prev = he1;

			he2->next = he3;
			he3->prev = he2;

			he3->next = he1;
			he1->prev = he3;

			// Set face
			myFace* face = new myFace();
			face->adjacent_halfedge = he1;
			he1->adjacent_face = face;
			he2->adjacent_face = face;
			he3->adjacent_face = face;

			// Set twins
			if (last_boundary_edge != nullptr)
			{
				he3->twin = last_boundary_edge;
				last_boundary_edge->twin = he3;
			}
			last_boundary_edge = he1;

			// Store edges and face
			faces.push_back(face);
			halfedges.push_back(he1);
			halfedges.push_back(he2);
			halfedges.push_back(he3);
		}
	}
	vertices.insert(vertices.end(), vertices1.begin(), vertices1.end());
	vertices.insert(vertices.end(), vertices2.begin(), vertices2.end());
}

bool myMesh::readFile(std::string filename)
{
	string s, t, u;
	vector<int> faceids;
	myHalfedge **hedges;

	ifstream fin(filename);
	if (!fin.is_open()) 
	{
		std::cout << "Unable to open file!\n";
		return false;
	}
	name = filename;

	map<pair<int, int>, myHalfedge *> halfEdge_map;
	map<pair<int, int>, myHalfedge *>::iterator it;
	
	while (getline(fin, s))
	{
		stringstream myline(s);
		myline >> t;
		if (t == "g") {}
		else if (t == "v")
		{
			float x, y, z;
			myline >> x >> y >> z;

			myVertex *v = new myVertex;
			v->point = new myPoint3D(x, y, z);
			v->index = vertices.size();
			vertices.push_back(v);
			//std::cout << "v " << x << " " << y << " " << z << " OK" << endl;
		}
		else if (t == "mtllib") {}
		else if (t == "usemtl") {}
		else if (t == "s") {}
		else if (t == "f")
		{
			//std::cout << "f"; 

			int i;
			int prev_i = -1;
			int first_i = -1;

			myHalfedge* prev_he = nullptr;
			myHalfedge* first_he = nullptr;

			myFace* face = new myFace;
			face->index = faces.size();
			faces.push_back(face);

			while (myline >> u)
			{
				if (prev_i == -1)
				{
					prev_i = atoi((u.substr(0, u.find("/"))).c_str()) - 1;
					first_i = prev_i;
					//std::cout << " " << first_i << " ";
					continue;
				}

				i = atoi((u.substr(0, u.find("/"))).c_str()) - 1;
				//std::cout << i << " ";

				// Create Half edges
				myHalfedge* edge = new myHalfedge;
				bool edge_exists = !halfEdge_map.insert(make_pair(make_pair(prev_i, i), edge)).second;

				// Set source
				edge->source = vertices[prev_i];
				edge->source->originof = edge;

				// Set prev
				if (prev_he != nullptr)
				{
					prev_he->next = edge;
					edge->prev = prev_he;
				}
				
				// Set face
				edge->adjacent_face = face;

				// Set first & prev data
				if (first_he == nullptr)
				{
					face->adjacent_halfedge = edge;
					first_he = edge;
				}
				prev_he = edge;
				prev_i = i;
			}

			// Connect last and first edges

			// Create Half edges
			myHalfedge* edge = new myHalfedge;
			bool edge_exists = !halfEdge_map.insert(make_pair(make_pair(i, first_i), edge)).second;

			// Set prev
			prev_he->next = edge;
			edge->prev = prev_he;

			// Set next
			edge->next = first_he;
			first_he->prev = edge;

			// Set source
			edge->source = vertices[i];
			edge->source->originof = edge;

			// Set face
			edge->adjacent_face = face;

			//std::cout << endl;
		}
	}

	for (const auto& kv : halfEdge_map) 
	{
		std::pair<int, int> key = kv.first;
		myHalfedge* value = kv.second;

		// Set all twins
		std::pair<int, int> twin_key = make_pair(key.second, key.first);
		if (halfEdge_map.find(twin_key) != halfEdge_map.end())
		{
			halfEdge_map[twin_key]->twin = value;
			value->twin = halfEdge_map[twin_key];
		}
		halfedges.push_back(value);
	}

	checkMesh();
	normalize();

	return true;
}

void myMesh::computeNormals()
{
	// Faces
	for (int i = 0; i < faces.size(); i++)
	{
		faces[i]->computeNormal();
	}

	// Vertices
	for (int i = 0; i < vertices.size(); i++)
	{
		vertices[i]->computeNormal();
	}
}

void myMesh::normalize()
{
	if (vertices.size() < 1) return;

	int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i]->point->X < vertices[tmpxmin]->point->X) tmpxmin = i;
		if (vertices[i]->point->X > vertices[tmpxmax]->point->X) tmpxmax = i;

		if (vertices[i]->point->Y < vertices[tmpymin]->point->Y) tmpymin = i;
		if (vertices[i]->point->Y > vertices[tmpymax]->point->Y) tmpymax = i;

		if (vertices[i]->point->Z < vertices[tmpzmin]->point->Z) tmpzmin = i;
		if (vertices[i]->point->Z > vertices[tmpzmax]->point->Z) tmpzmax = i;
	}

	double xmin = vertices[tmpxmin]->point->X, xmax = vertices[tmpxmax]->point->X,
		ymin = vertices[tmpymin]->point->Y, ymax = vertices[tmpymax]->point->Y,
		zmin = vertices[tmpzmin]->point->Z, zmax = vertices[tmpzmax]->point->Z;

	double scale = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
	scale = scale > (zmax - zmin) ? scale : (zmax - zmin);

	for (unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i]->point->X -= (xmax + xmin) / 2;
		vertices[i]->point->Y -= (ymax + ymin) / 2;
		vertices[i]->point->Z -= (zmax + zmin) / 2;

		vertices[i]->point->X /= scale;
		vertices[i]->point->Y /= scale;
		vertices[i]->point->Z /= scale;
	}
}

void myMesh::splitFaceTRIS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}

void myMesh::splitEdge(myHalfedge *e1, myPoint3D *p)
{
	/**** TODO ****/
}

void myMesh::splitFaceQUADS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}

std::set<myVertex*> getFaceVertices(myFace* face)
{
	std::set<myVertex*> vertices;
	myHalfedge* he = face->adjacent_halfedge;
	do
	{
		vertices.insert(he->source);
		he = he->next;
	} while (he != face->adjacent_halfedge);
	return vertices;
}

std::set<myFace*> getNeighborFaces(myVertex* v)
{
	std::set<myFace*> faces;

	myHalfedge* current_he = v->originof;
	do
	{
		faces.insert(current_he->adjacent_face);
		current_he = current_he->twin->next;
	} 
	while (current_he != v->originof);

	return faces;
}

std::set<myHalfedge*> getNeighborEdges(myVertex* v)
{
	std::set<myHalfedge*> edges;

	myHalfedge* current_he = v->originof;
	do
	{
		edges.insert(current_he);
		current_he = current_he->twin->next;
	} while (current_he != v->originof);

	return edges;
}

myPoint3D getAveragePoint(std::set<myPoint3D> points)
{
	myPoint3D res = myPoint3D(0, 0, 0);
	for (myPoint3D p : points)
	{
		res += p;
	}
	res /= static_cast<double>(points.size());
	return res;
}

myPoint3D getAveragePoint(std::set<myVertex*> vertices)
{
	std::set<myPoint3D> points;
	for (myVertex* v : vertices)
	{
		points.insert(*(v->point));
	}
	return getAveragePoint(points);
}

myPoint3D getBarycenter(myPoint3D a, myPoint3D b, myPoint3D c, double wa, double wb, double wc, int n)
{
	return (a * wa + b * wb + c * wc) / n;
}

/// <summary>
/// https://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface
/// </summary>
void myMesh::subdivisionCatmullClark()
{
	//TODO: Check source, originof, next, prev, twin, adjacent_face, adjacent_he

	std::map<myFace*, myVertex*> face_vertex_map;
	std::map<myHalfedge*, myVertex*> edge_vertex_map;
	std::map <std::pair<myVertex*, myVertex*>, myHalfedge*> he_map;

	std::vector<myHalfedge*> new_edges;
	std::vector<myFace*> new_faces;
	std::vector<myVertex*> new_vertices;

	// Generate face points
	for (myFace* face : faces)
	{
		std::set<myVertex*> face_vertices = getFaceVertices(face);
		myVertex* face_vertex = new myVertex();
		face_vertex->point = new myPoint3D((getAveragePoint(face_vertices)));
		face_vertex_map.insert(std::pair<myFace*, myVertex*>(face, face_vertex));
		new_vertices.push_back(face_vertex);
	}

	// Generate edge points
	for (myHalfedge* he : halfedges)
	{
		if (he < he->twin) // Avoid duplicate vertices
		{
			std::set<myVertex*> vertices_used;
			vertices_used.insert(face_vertex_map.at(he->adjacent_face));
			vertices_used.insert(face_vertex_map.at(he->twin->adjacent_face));
			vertices_used.insert(he->source);
			vertices_used.insert(he->twin->source);

			myVertex* edge_vertex = new myVertex();
			edge_vertex->point = new myPoint3D(getAveragePoint(vertices_used));
			edge_vertex_map.insert(std::pair<myHalfedge*, myVertex*>(he, edge_vertex));
			edge_vertex_map.insert(std::pair<myHalfedge*, myVertex*>(he->twin, edge_vertex));
			new_vertices.push_back(edge_vertex);
		}
	}

	// Average original points
	for (myVertex* v : vertices)
	{
		std::set<myFace*> neighbor_faces = getNeighborFaces(v);
		std::set<myHalfedge*> neighbor_edges = getNeighborEdges(v);

		// Average of face vertices
		std::set<myVertex*> face_vertices;
		for (myFace* f : neighbor_faces)
		{
			face_vertices.insert(face_vertex_map.at(f));
		}
		myPoint3D avg_faces = getAveragePoint(face_vertices);

		// Average of edge vertices
		std::set<myPoint3D> edge_middle_points;
		for (myHalfedge* he : neighbor_edges)
		{
			std::set<myPoint3D> he_points;
			he_points.insert(*(he->source->point));
			he_points.insert(*(he->twin->source->point));
			edge_middle_points.insert(getAveragePoint(he_points));
		}
		myPoint3D avg_edges = getAveragePoint(edge_middle_points);

		// Average point
		int n = neighbor_faces.size();
		myPoint3D barycenter = getBarycenter(avg_faces, avg_edges, *(v->point), 1, 2, n - 3, n);
		
		if (barycenter.dist(myPoint3D(0, 0, 0)) > 1.0)
		{
			std::cout << "breakpoint" << std::endl;
		}

		v->point->copyValuesFrom(barycenter);		

	}

	// Create edgeVertex <-> faceVertex half edges
	for (std::pair<myFace*, myVertex*> f_v_pair : face_vertex_map)
	{
		myFace* face = f_v_pair.first;
		myVertex* face_vertex = f_v_pair.second;

		std::set<myVertex*> edge_vertices;
		for (std::pair<myHalfedge*, myVertex*> he_v_pair : edge_vertex_map)
		{
			myHalfedge* he = he_v_pair.first;
			myVertex* edge_vertex = he_v_pair.second;

			if (he->adjacent_face == face)
			{
				edge_vertices.insert(edge_vertex);
			}
		}
		
		for (myVertex* edge_vertex : edge_vertices)
		{
			myHalfedge* he_from_face_vertex = new myHalfedge();
			he_from_face_vertex->source = face_vertex;
			face_vertex->originof = he_from_face_vertex;

			myHalfedge* he_to_face_vertex = new myHalfedge();
			he_to_face_vertex->source = edge_vertex;
			edge_vertex->originof = he_to_face_vertex;

			he_from_face_vertex->twin = he_to_face_vertex;
			he_to_face_vertex->twin = he_from_face_vertex;

			// Fill map
			auto mapping1 = make_pair(make_pair(face_vertex, edge_vertex), he_from_face_vertex);
			auto mapping2 = make_pair(make_pair(edge_vertex, face_vertex), he_to_face_vertex);

			he_map.insert(mapping1);
			he_map.insert(mapping2);

			new_edges.push_back(he_from_face_vertex);
			new_edges.push_back(he_to_face_vertex);
		}
	}
	
	// Create originalVertex <-> edgeVertex half edges
	for (myVertex* original_v : vertices)
	{
		std::set<myVertex*> edge_vertices;
		for (myHalfedge* he : halfedges)
		{
			if (he->source == original_v)
			{
				edge_vertices.insert(edge_vertex_map.at(he));
			}
		}
		
		for (myVertex* edge_vertex : edge_vertices)
		{
			myHalfedge* he_to_edge = new myHalfedge();
			he_to_edge->source = original_v;
			original_v->originof = he_to_edge;

			myHalfedge* he_from_edge = new myHalfedge();
			he_from_edge->source = edge_vertex;
			edge_vertex->originof = he_from_edge;

			he_to_edge->twin = he_from_edge;
			he_from_edge->twin = he_to_edge;

			// Fill map
			auto mapping1 = make_pair(make_pair(original_v, edge_vertex), he_to_edge);
			auto mapping2 = make_pair(make_pair(edge_vertex, original_v), he_from_edge);

			he_map.insert(mapping1);
			he_map.insert(mapping2);

			new_edges.push_back(he_from_edge);
			new_edges.push_back(he_to_edge);
		}
	}

	// Connect all half edges
	for (myFace* face : faces)
	{
		myHalfedge* current_he = face->adjacent_halfedge;
		do
		{
			myVertex* face_vertex = face_vertex_map.at(face);
			myVertex* original_vertex = current_he->source;
			myVertex* edge_vertex = edge_vertex_map.at(current_he);
			myVertex* prev_edge_vertex = edge_vertex_map.at(current_he->prev);

			// Connect edges
			he_map.at(make_pair(original_vertex, edge_vertex))->next = he_map.at(make_pair(edge_vertex, face_vertex));
			he_map.at(make_pair(edge_vertex, face_vertex))->prev = he_map.at(make_pair(original_vertex, edge_vertex));

			he_map.at(make_pair(edge_vertex, face_vertex))->next = he_map.at(make_pair(face_vertex, prev_edge_vertex));
			he_map.at(make_pair(face_vertex, prev_edge_vertex))->prev = he_map.at(make_pair(edge_vertex, face_vertex));

			he_map.at(make_pair(face_vertex, prev_edge_vertex))->next = he_map.at(make_pair(prev_edge_vertex, original_vertex));
			he_map.at(make_pair(prev_edge_vertex, original_vertex))->prev = he_map.at(make_pair(face_vertex, prev_edge_vertex));

			he_map.at(make_pair(prev_edge_vertex, original_vertex))->next = he_map.at(make_pair(original_vertex, edge_vertex));
			he_map.at(make_pair(original_vertex, edge_vertex))->prev = he_map.at(make_pair(prev_edge_vertex, original_vertex));

			// Create face
			myFace* new_face = new myFace();
			he_map.at(make_pair(original_vertex, edge_vertex))->adjacent_face = new_face;
			he_map.at(make_pair(edge_vertex, face_vertex))->adjacent_face = new_face;
			he_map.at(make_pair(face_vertex, prev_edge_vertex))->adjacent_face = new_face;
			he_map.at(make_pair(prev_edge_vertex, original_vertex))->adjacent_face = new_face;
			new_face->adjacent_halfedge = he_map.at(make_pair(original_vertex, edge_vertex));
			new_faces.push_back(new_face);

			current_he = current_he->next;
		}
		while (current_he != face->adjacent_halfedge);
	}

	// Remove old halfedges and faces
	for (myFace* f : faces)
	{
		freeFace(f);
	}
	faces.clear();

	for (myHalfedge* he : halfedges)
	{
		freeHalfEdge(he);
	}
	halfedges.clear();
	
	vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
	faces = new_faces;
	halfedges = new_edges;
	checkMesh();
}

bool isTriangle(myFace* f)
{
	// Check if face is a triangle
	if (f == nullptr || f->adjacent_halfedge == nullptr || f->adjacent_halfedge->next == nullptr || f->adjacent_halfedge->next->next == nullptr || f->adjacent_halfedge->next->next->next == nullptr) return false;
	return f->adjacent_halfedge->next->next->next == f->adjacent_halfedge;
}

void myMesh::triangulate()
{
	std::vector<myFace*> temp_faces = faces;
	int fails = 0;
	for (myFace* f : temp_faces)
	{
		if (!triangulate(f))
		{
			fails++;
		}
	}
	if (fails > 0)
	{
		std::cerr << "Failed to triangulate " << fails << " face(s)!" << std::endl;
	}
	checkMesh();
}

myHalfedge* myMesh::getShortestValidEdgeNonTriangle()
{
	double min = INFINITY;
	myHalfedge* min_he = nullptr;
	int min_n_edges = INFINITY;
	for (myHalfedge* he : halfedges)
	{
		// Check if face is a triangle
		int n_edges = 0;
		bool is_valid = true;
		myHalfedge* current_he = he;
		do 
		{
			if (current_he->twin == nullptr)
			{
				is_valid = false;
				break;
			}
			current_he = current_he->next;
			n_edges++;
		}
		while (current_he != he);

		// Check if twin face is a triangle
		int n_twin_edges = 0;
		current_he = he->twin;
		do
		{
			if (current_he->twin == nullptr)
			{
				is_valid = false;
				break;
			}
			current_he = current_he->next;
			n_twin_edges++;
		} 
		while (current_he != he->twin);

		if (!is_valid || n_edges <= 3 || n_twin_edges <= 3) continue;

		myPoint3D* p1 = he->source->point;
		myPoint3D* p2 = he->twin->source->point;

		double dist = p1->dist(*p2);
		if (dist < min)
		{
			min = dist;
			min_he = he;
			min_n_edges = n_edges;
		}
	}
	return min_he;
}

myHalfedge* myMesh::getShortestValidEdge()
{
	double min = INFINITY;
	myHalfedge* min_he = nullptr;
	int min_n_edges = INFINITY;
	for (myHalfedge* he : halfedges)
	{
		int n_edges = 0;
		bool is_valid = true;
		myHalfedge* current_he = he;
		do
		{
			if (he->twin == nullptr)
			{
				is_valid = false;
				break;
			}
			he = he->next;
			n_edges++;
		} while (he != current_he);

		if (!is_valid) continue;

		myPoint3D* p1 = he->source->point;
		myPoint3D* p2 = he->twin->source->point;

		double dist = p1->dist(*p2);
		if (dist < min)
		{
			min = dist;
			min_he = he;
			min_n_edges = n_edges;
		}
	}
	return min_he;
}

void myMesh::simplify()
{
	myHalfedge* shortest_edge = getShortestValidEdge();
	if (shortest_edge == nullptr || faces.size() < 3)
	{
		std::cerr << "No polygon left for simplification" << std::endl;
	}
	else if (!unifyEdge(shortest_edge, COLLAPSE_AVERAGE))
	{
		std::cerr << "Abort simplification because he->source is equal to he->next->source" << std::endl;
	}
}

void myMesh::resetOriginOf()
{
	for (myHalfedge* he : halfedges)
	{
		he->source->originof = he;
	}
}

void myMesh::freeVertex(myVertex* v)
{
	for (myHalfedge* he : halfedges)
	{
		if (he->source == v)
		{
			he->source = nullptr;
		}
	}
	delete v;
}

void myMesh::freePoint(myPoint3D* p)
{

	for (myVertex* v : vertices)
	{
		if (v->point == p)
		{
			v->point = nullptr;
		}
	}
	delete p;
}

void myMesh::freeHalfEdge(myHalfedge* he)
{
	for (myVertex* v : vertices)
	{
		if (v->originof == he)
		{
			v->originof = nullptr;
		}
	}

	for (myHalfedge* m_he : halfedges)
	{
		if (m_he->next == he)
		{
			m_he->next = nullptr;
		}
		if (m_he->prev == he)
		{
			m_he->prev = nullptr;
		}
		if (m_he->twin == he)
		{
			m_he->twin = nullptr;
		}
	}

	for (myFace* f : faces)
	{
		if (f->adjacent_halfedge == he)
		{
			f->adjacent_halfedge = nullptr;
		}
	}
	delete he;
}

void myMesh::freeFace(myFace* f)
{
	if (f->adjacent_halfedge != nullptr)
	{
		f->adjacent_halfedge->adjacent_face = nullptr;
	}
	for (myHalfedge* he : halfedges)
	{
		if (he->adjacent_face == f)
		{
			he->adjacent_face = nullptr;
		}
	}
	delete f;
}

bool myMesh::unifyEdge(myHalfedge* he, CollapseMode mode)
{
	myHalfedge* he_to_delete = he;
	myHalfedge* twin_to_delete = he->twin;

	myVertex* v_to_delete = he->source;
	myVertex* v_target = he->next->source;
	v_target->originof = he->next;
	

	bool face_is_triangle = he->next->next->next == he;
	bool twin_is_triangle = he->twin->next->next->next == he->twin;

	// Use if not willing to simplify triangles
	//if (face_is_triangle) return false;
	//if (twin_is_triangle) return false;

	if (v_target == v_to_delete) return false;

	myHalfedge* twin_A1 = he->next->twin;
	myHalfedge* twin_B2 = he->twin->prev->twin;

	myHalfedge* twin_1A = he->prev->twin;
	myHalfedge* twin_2B = he->twin->next->twin;

	he->adjacent_face->adjacent_halfedge = he->next;
	he->twin->adjacent_face->adjacent_halfedge = he->twin->next;

	// Set sources and originofs
	myHalfedge* current_he = he;
	do
	{
		current_he = current_he->twin->next;
		current_he->source = v_target;
		current_he->source->originof = current_he;
	}
	while (current_he != he);

	twin_A1->source->originof = twin_A1;
	twin_B2->source->originof = twin_B2;

	he->prev->next = he->next;
	he->next->prev = he->prev;

	he->twin->prev->next = he->twin->next;
	he->twin->next->prev = he->twin->prev;

	// Move vertex
	if (mode == COLLAPSE_START)
	{
		v_target->point->X = v_to_delete->point->X;
		v_target->point->Y = v_to_delete->point->Y;
		v_target->point->Z = v_to_delete->point->Z;
	}
	else if (mode == COLLAPSE_AVERAGE)
	{
		v_target->point->X = (v_target->point->X + v_to_delete->point->X) / 2;
		v_target->point->Y = (v_target->point->Y + v_to_delete->point->Y) / 2;
		v_target->point->Z = (v_target->point->Z + v_to_delete->point->Z) / 2;
	}
	else if (mode == COLLAPSE_END)
	{

	}
	
	if (face_is_triangle)
	{
		twin_A1->twin = twin_1A;
		twin_1A->twin = twin_A1;

		faces.erase(std::remove(faces.begin(), faces.end(), he_to_delete->adjacent_face), faces.end());
		freeFace(he_to_delete->adjacent_face);

		halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), he_to_delete->next->next), halfedges.end());
		freeHalfEdge(he_to_delete->next->next);

		halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), he_to_delete->next), halfedges.end());
		freeHalfEdge(he_to_delete->next);
	}
	if (twin_is_triangle)
	{
		twin_B2->twin = twin_2B;
		twin_2B->twin = twin_B2;

		faces.erase(std::remove(faces.begin(), faces.end(), twin_to_delete->adjacent_face), faces.end());
		freeFace(twin_to_delete->adjacent_face);

		halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), twin_to_delete->next->next), halfedges.end());
		freeHalfEdge(twin_to_delete->next->next);

		halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), twin_to_delete->next), halfedges.end());
		freeHalfEdge(twin_to_delete->next);
	}

	halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), he_to_delete), halfedges.end());
	freeHalfEdge(he_to_delete);

	halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), twin_to_delete), halfedges.end());
	freeHalfEdge(twin_to_delete);

	vertices.erase(std::remove(vertices.begin(), vertices.end(), v_to_delete), vertices.end());
	freeVertex(v_to_delete);

	// Set originofs 
	for (myHalfedge* he_to_check : halfedges)
	{
		if (!he_to_check->source) continue;
		he_to_check->source->originof = he_to_check;
	}

	return true;
}

bool myMesh::triangulate(myFace *f)
{
	// Check if face is a triangle
	if (isTriangle(f))
	{
		return false;
	}

	// Calculate the center of the face
	myPoint3D center(0, 0, 0);
	int count = 0;
	myHalfedge* he = f->adjacent_halfedge;
	do
	{
		center += *(he->source->point);
		he = he->next;
		count++;
	} 
	while (he != f->adjacent_halfedge);
	center /= count;

	// Create a new vertex at the center
	myVertex* center_vertex = new myVertex;
	center_vertex->point = new myPoint3D(center);
	vertices.push_back(center_vertex);

	// Store the original half-edges
	std::vector<myHalfedge*> original_halfedges;
	he = f->adjacent_halfedge;
	do
	{
		original_halfedges.push_back(he);
		he = he->next;
	} 
	while (he != f->adjacent_halfedge);

	// Store new half-edges
	std::vector<myHalfedge*> hes_to_center;
	std::vector<myHalfedge*> hes_from_center;

	// Create new faces and half-edges
	for (myHalfedge* he : original_halfedges)
	{
		myHalfedge* he_from_center = new myHalfedge;
		myHalfedge* he_to_center = new myHalfedge;

		if (center_vertex->originof == nullptr)
		{
			center_vertex->originof = he_from_center;
		}

		myFace* new_face = new myFace;
		new_face->adjacent_halfedge = he_from_center;
		he_from_center->adjacent_face = new_face;
		he_to_center->adjacent_face = new_face;

		//TODO: Fix non sense (inverted hes)
		he_to_center->source = he->next->source;
		he_to_center->source->originof = he_to_center;
		he_to_center->next = he_from_center;
		he_to_center->prev = he;

		he_from_center->source = center_vertex;
		he_from_center->source->originof = he_from_center;
		he_from_center->next = he;
		he_from_center->prev = he_to_center;
		he->next = he_to_center;
		he->prev = he_from_center;
		he->adjacent_face = new_face;

		faces.push_back(new_face);
		halfedges.push_back(he_to_center);
		halfedges.push_back(he_from_center);
		hes_from_center.push_back(he_from_center);
		hes_to_center.push_back(he_to_center);
	}

	// Update twins
	for (size_t i = 0; i < hes_to_center.size(); ++i)
	{
		hes_to_center[i]->twin = hes_from_center[(i + 1) % hes_from_center.size()];
		hes_from_center[i]->twin = hes_to_center[(i - 1 + hes_to_center.size()) % hes_to_center.size()];
	}

	// Remove the original face
	faces.erase(std::remove(faces.begin(), faces.end(), f), faces.end());
	delete f;

	int i = 0;
	for (myHalfedge* he : halfedges) he->index = i++;
	return true;
}

