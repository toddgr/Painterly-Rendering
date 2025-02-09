/*

Data structures for learnply

Eugene Zhang 2005

*/

#ifndef __LEARNPLY_H__
#define __LEARNPLY_H__


#include "ply.h"
#include "icVector.H"

const double EPS = 1.0e-6;
const double PI = 3.1415926535898;

/* forward declarations */
class Quad;
class Edge;

class Vertex {
public:
	double x, y, z;			/*coordinates*/

	double vx, vy, vz;		/*vector field*/
	double scalar = 0;		/*scalar field*/

	double R = 0, G = 0, B = 0;		/*color*/

	int index;

	int nquads;
	Quad** quads;
	int max_quads;

	int nedges;
	Edge** edges;
	int max_edges;

	icVector3 normal;
	void* other_props = NULL;

	//Vertex& operator= (const Vertex& b) {
	//	this->x = b.x;
	//	this->y = b.y;
	//	this->z = b.z;
	//	this->R = b.R;
	//	this->G = b.G;
	//	this->B = b.B;

	//	return *this;
	//}
public:
	Vertex(double xx, double yy, double zz) { 
		x = xx; 
		y = yy; 
		z = zz; 
	}

	Vertex(double xx, double yy, double zz, double rr, double gg, double bb) {
		x = xx;
		y = yy;
		z = zz;
		R = rr;
		G = gg;
		B = bb;
	}
};

class Edge {
public:
	int index;

	Vertex* verts[2];

	int nquads;
	Quad** quads;

	double length;
};

class Quad {
public:
	int index;

	Vertex* verts[4];
	Edge* edges[4];

	float area;

	icVector3 normal;
	void* other_props;

	icVector2* singularity;
};


class Polyhedron {
public:

	Quad** qlist;		/* list of quads */
	int nquads;
	int max_quads;

	Vertex** vlist;    /* list of vertices */
	int nverts;
	int max_verts;

	Edge** elist;      /* list of edges */
	int nedges;
	int max_edges;

	icVector3 center;
	double radius;
	double area;

	int selected_quad;
	int selected_vertex;
	unsigned char orientation;  // 0=ccw, 1=cw

	PlyOtherProp* vert_other, * face_other;

	/*constructors*/
	Polyhedron();
	Polyhedron(FILE*);

	/*initialization functions*/
	void create_pointers();
	void average_normals();
	void create_edge(Vertex*, Vertex*);
	void create_edges();
	int face_to_vertex_ref(Quad*, Vertex*);
	void order_vertex_to_quad_ptrs(Vertex*);
	void vertex_to_quad_ptrs();
	void vertex_to_edge_ptrs();
	void calc_bounding_sphere();
	void calc_face_normals_and_area();
	void calc_edge_length();
	Vertex* other_vert(Vertex* vert, Edge* edge);
	Edge* find_edge(Vertex* v1, Vertex* v2);
	Quad* find_quad(double x, double y);
	double smallest_x(Quad* temp);
	double largest_x(Quad* temp);
	double smallest_y(Quad* temp);
	double largest_y(Quad* temp);

	/*utilties*/
	Quad* find_common_edge(Quad*, Vertex*, Vertex*);
	Quad* other_quad(Edge*, Quad*);

	/*feel free to add more to help youself*/
	void write_info();
	void write_file(FILE*);


	/*initialization and finalization*/
	void initialize();
	void finalize();
};

#endif /* __LEARNPLY_H__ */

