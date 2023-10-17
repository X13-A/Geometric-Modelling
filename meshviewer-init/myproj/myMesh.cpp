#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myvector3d.h"

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
	vector<myHalfedge *>::iterator it;

	int twin_errors = 0;
	int next_errors = 0;
	int prev_errors = 0;
	int source_errors = 0;
	int adjacent_face_errors = 0;

	for (it = halfedges.begin(); it != halfedges.end(); it++)
	{
		if ((*it)->twin == NULL) twin_errors++;
		if ((*it)->next == NULL) next_errors++;
		if ((*it)->prev == NULL) prev_errors++;
		if ((*it)->source == NULL) source_errors++;
		if ((*it)->adjacent_face == NULL) adjacent_face_errors++;
	}
	if (twin_errors != 0)
	{
		std::cout << "Error! Not all edges have their twin (" << twin_errors << ")!\n";
	}
	if (next_errors != 0)
	{
		std::cout << "Error! Not all edges have their next (" << next_errors << ")!\n";
	}
	if (prev_errors != 0)
	{
		std::cout << "Error! Not all edges have their prev (" << prev_errors << ")!\n";
	}
	if (source_errors != 0)
	{
		std::cout << "Error! Not all edges have their source (" << source_errors << ")!\n";
	}
	if (adjacent_face_errors != 0)
	{
		std::cout << "Error! Not all edges have their adjacent_face (" << adjacent_face_errors << ")!\n";
	}
	if (twin_errors + next_errors + prev_errors + source_errors + adjacent_face_errors == 0) std::cout << "Each edge is set!\n";
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

void myMesh::subdivisionCatmullClark()
{
	/**** TODO ****/
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

void myMesh::simplify()
{
	triangulate();
	if (!unifyEdge(halfedges[5], COLLAPSE_AVERAGE))
	{
		std::cerr << "Failed simplification" << std::endl;
	}
	std::cerr << vertices.size() << " Vertices" << std::endl;
	std::cerr << faces.size() << " Faces" << std::endl;
	std::cerr << halfedges.size() << " Half edges" << std::endl;

	resetOriginOf();
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
	if (v->originof)
	{
		v->originof->source = nullptr;
	}
	for (myHalfedge* he : halfedges)
	{
		if (he->source == v)
		{
			he->source = nullptr;
		}
	}

	vertices.erase(std::remove(vertices.begin(), vertices.end(), v), vertices.end());
	delete v;
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

	halfedges.erase(std::remove(halfedges.begin(), halfedges.end(), he), halfedges.end());
	delete he;
}

void myMesh::freeFace(myFace* f)
{
	if (f->adjacent_halfedge)
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
	faces.erase(std::remove(faces.begin(), faces.end(), f), faces.end());
	delete f;
}

bool myMesh::unifyEdge(myHalfedge* he, CollapseMode mode)
{
	// STEPS
	// - TODO: Check if unification is possible
	// - Mark vertex to delete & new vertex
	// - Mark half-edges to delete
	// - Set neighbors source to new vertex
	// - Set new twins to replace deleted half-edges
	// - Delete vertex, faces & half-edges

	if (he == nullptr) return false;
	if (he->twin == nullptr) return false;
	if (he->prev == nullptr) return false;
	if (he->next == nullptr) return false;

	myHalfedge* twin_A1 = he->prev->twin;
	myHalfedge* twin_B2 = he->twin->next->twin;

	myHalfedge* twin_1A = he->next->twin;
	myHalfedge* twin_2B = he->twin->prev->twin;

	if (!twin_A1 || !twin_1A || !twin_B2 || !twin_2B) return false;

	// Mark half-edges defining the triangles to delete
	std::vector<myHalfedge*> hes_to_delete =
	{
		// First triangle
		he, he->next, he->next->next, 
		// Second triangle
		he->twin, he->twin->next, he->twin->next->next
	};

	for (myHalfedge* current_he : hes_to_delete)
	{
		if (!current_he) return false;
	}

	myVertex* v_to_delete = he->source;
	myVertex* v_target = he->next->source;


	// Use new vertex
	myHalfedge* current_he = he;
	do
	{
		current_he = current_he->twin->next;
		current_he->source = v_target;
		current_he->source->originof = current_he;
	}
	while (current_he != he);

	// Set new twins
	twin_A1->twin = twin_1A;
	twin_1A->twin = twin_A1;

	twin_B2->twin = twin_2B;
	twin_2B->twin = twin_B2;

	// Set originof
	twin_A1->twin->source->originof = twin_A1->twin;
	twin_B2->source->originof = twin_B2;

	current_he = twin_A1;
	do
	{
		current_he->source = v_target;
		current_he->twin->source->originof = current_he->twin;
		current_he = current_he->twin->next;
	}
	while (current_he != twin_A1);
	v_target->originof = twin_A1;

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

	// Delete old data
	freeVertex(v_to_delete);
	std::vector<myFace*> faces_to_delete = {};

	for (myHalfedge* he_to_delete : hes_to_delete)
	{
		// Mark faces for deletion
		myFace* f = he_to_delete->adjacent_face;
		if (std::find(faces_to_delete.begin(), faces_to_delete.end(), f) == faces_to_delete.end())
		{
			faces_to_delete.push_back(f);
		}
		freeHalfEdge(he_to_delete);
	}

	for (myFace* f : faces_to_delete)
	{
		freeFace(f);
	}

	return true;
}
 
//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace *f)
{
	// Check if face is a triangle
	if (f == nullptr || f->adjacent_halfedge == nullptr || f->adjacent_halfedge->next == nullptr || f->adjacent_halfedge->next->next == nullptr || f->adjacent_halfedge->next->next->next == nullptr) return false;
	if (f->adjacent_halfedge->next->next->next == f->adjacent_halfedge)
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

