#include "VectorFieldTopology.h"
#include "Sobel.h"
#include <iostream>
#define M_PI 3.14159265358979323846
#define NPN 256 // for sobel-- need to change this later

constexpr auto EPSILON = 1.0e-5;
constexpr auto STEP = 0.005;
constexpr auto MIN_K = 0.05;

// min and max texture coords
int rmin = NPN;
int rmax = 0;
int cmin = 0;
int cmax = NPN; 

extern Polyhedron* poly;
extern std::vector<POLYLINE> polylines;
extern std::list<Singularity> singularities;
extern int win_width;
extern GLubyte patsvec[NPN][NPN][2];

icVector3 min, max;

bool isZero(double x) {
	if (x == 0.0) {
		return true;
	}
	else {
		return false;
	}
}


bool sinp2Boundary(icVector3& currPos, const icVector3& min, const icVector3& max) {
	bool hitBoundary = false;
	if (currPos.x < min.x) {
		hitBoundary = true;
		currPos.x = min.x;
	}
	if (currPos.y < min.y) {
		hitBoundary = true;
		currPos.y = min.y;
	}
	if (currPos.x > max.x) {
		hitBoundary = true;
		currPos.x = max.x;
	}
	if (currPos.y > max.y) {
		hitBoundary = true;
		currPos.y = max.y;
	}
	return hitBoundary;
}

// Streamline tracing
// Traces the streamline into the next quad
void streamlineTrace(Quad*& nextQuad, Quad* currQuad, icVector3 currPos, icVector3 currVec, double t,
	const icVector3& min, const icVector3& max) {

	bool insideQuad = false;
	while (!insideQuad) {

		// Is outside the field
		if (currPos.x < min.x || currPos.x > max.x ||
			currPos.y < min.y || currPos.y > max.y) {
			nextQuad = nullptr;
			return;
		}

		double t_ = INFINITY;
		Quad* nextQuad_ = nullptr;
		for (int e_i = 0; e_i < 4; e_i++) { // For each of the quad edges
			Edge* edge = currQuad->edges[e_i];
			Vertex* v0 = edge->verts[0];
			Vertex* v1 = edge->verts[1];
			double t_temp;
			if (std::abs(v0->x - v1->x) < EPSILON) {
				t_temp = (v0->x - currPos.x) / currVec.x;
			}
			else {
				t_temp = (v0->y - currPos.y) / currVec.y;
			}
			if (t_temp > 0 && t_temp < t_) {
				// Get next quad
				t_ = t_temp;
				if (edge->quads[0] != currQuad && edge->quads[1] == currQuad) {
					nextQuad_ = edge->quads[0];
				}
				else if (edge->quads[0] == currQuad && edge->quads[1] != currQuad) {
					nextQuad_ = edge->quads[1];
				}
			}
		}
		if (nextQuad_ == nullptr) {
			/*t_ = t / 1000;
			t = t - t_;
			currPos = currPos + currVec * t_;
			currQuad = findQuad(currPos);
			*/
			currPos = currPos + currVec * t;
			nextQuad = findQuad(currPos);
			return;
		}
		else {
			if (t_ >= t) {
				insideQuad = true;
			}
			else {
				currQuad = nextQuad_;
				t = t - t_;
				currPos = currPos + currVec * t_;
			}
		}
	}
	nextQuad = currQuad;
}

// Get streamline (either forward or backward) from a given vertex
void streamlineFB(POLYLINE& line, const icVector3& seed, const double& step, bool forward) {  // Verified
	line.m_vertices.push_back(seed);	// Push the initial vertex to the line
	Quad* quad = findQuad(seed);		// Find the quad that the initial vertex is in
	icVector3 currPos = seed;			// Advance to the next
	double coef = 1.0;
	if (!forward) {
		coef = -1.0;
	}
	while (quad != nullptr) {
		icVector3 currVec = getVector(quad, currPos);	// Get the Vector
		//v:=0											// Could this be replaced with the edge field?
		if (currVec.length() < EPSILON) {				// Control the length of the vector
			break;
		}
		//first euler
		icVector3 nextPos = currPos + step * currVec * coef;
		// Is the next position in the boundary?
		if (sinp2Boundary(nextPos, min, max)) {
			line.m_vertices.push_back(nextPos);	// If yes add to the line
			break;
		}
		//get the new quad
		Quad* nextQuad = nullptr;
		// Trace the streamline into the next quad
		streamlineTrace(nextQuad, quad, currPos, currVec * coef, step, min, max);
		//update quad, pos, vector
		quad = nextQuad;
		currPos = nextPos;
		line.m_vertices.push_back(currPos);
	}
}


void streamline(POLYLINE& line, const icVector3& seed, const double& step) {  // verified
	streamlineFB(line, seed, step);				// Create streamline forward
	POLYLINE line_back;
	streamlineFB(line_back, seed, step, false);	// Create streamline backward
	line.merge(line_back);						// Merge the two together
}


void drawstreamlines() {
	POLYLINE line;
	findMinMaxField(min, max);			// Find the minimum and maximum coordinates
	for (int i = -1; i < 1; i++) {	// Display streamlines
		line.m_vertices.clear();
		streamline(line, icVector3(i, i, 0), 0.1);	// d2 was 0.001 but was taking too long to render
		line.m_rgb = icVector3(1.0, 0.0, 0.0);			// Streamlines are white for now
		polylines.push_back(line);						// Add line to polylines
		printf("streamline drawn\n");
	}
}

// Find minimum and maximum coordinate
void findMinMaxField(icVector3& min, icVector3& max) {
	min.x = poly->vlist[0]->x;
	min.y = poly->vlist[0]->y;
	min.z = poly->vlist[0]->z;
	max = min;
	for (int i = 1; i < poly->nverts; i++) {
		if (min.x > poly->vlist[i]->x) {
			min.x = poly->vlist[i]->x;
		}
		if (min.y > poly->vlist[i]->y) {
			min.y = poly->vlist[i]->y;
		}
		if (min.z > poly->vlist[i]->z) {
			min.z = poly->vlist[i]->z;
		}

		if (max.x < poly->vlist[i]->x) {
			max.x = poly->vlist[i]->x;
		}
		if (max.y < poly->vlist[i]->y) {
			max.y = poly->vlist[i]->y;
		}
		if (max.z < poly->vlist[i]->z) {
			min.z = poly->vlist[i]->z;
		}
		
	}
	std::cout << "min: {" << min.x << ", " << min.y << ", " << min.z << "}" << std::endl;
	std::cout << "max: {" << max.x << ", " << max.y << ", " << max.z << "}" << std::endl;
}

Quad* findQuad(const icVector3& v) {  // verified
	for (int i = 0; i < poly->nquads; i++) {
		Quad* qtemp = poly->qlist[i];
		if (insideQuad(qtemp, v))
			return qtemp;
	}
	return nullptr;
}

// Inside quad
bool insideQuad(const Quad* q, const icVector3& p) {  // verified
	double v0x = q->verts[2]->x;
	double v0y = q->verts[2]->y;
	double v2x = q->verts[0]->x;
	double v2y = q->verts[0]->y;

	if (p.x >= v0x && p.x <= v2x && p.y >= v0y && p.y <= v2y) {
		return true;
	}
	else {
		return false;
	}
}


// Get the vector from vector field by *bilinear interpolation*
// Could we alter this to instead get the gradient vector from the edge field?
icVector3 getVector(Quad* q, const icVector3& p) {
	// min.x, max.x,
	// min.y, max.y
	// rmin, rmax
	// cmin, cmax
	double x0, x1, x2, x3,
		x0p, x1p, x2p, x3p,
		y0, y1, y2, y3,
		y0p, y1p, y2p, y3p,
		vc0,
		vr0,
		vx0, vx1, vx2, vx3,
		vy0, vy1, vy2, vy3;
	double vz = 0.; // Need to think about what to have for z vector
	int r0, r1, r2, r3,
		c0, c1, c2, c3;

	// Get the vertices from the quad space
	x0 = q->verts[0]->x;
	y0 = q->verts[0]->y;

	// Get the corresponding texels in the texture space
	c0 = ((cmax - cmin) / (max.x - min.x)) * x0 + (((cmin * max.x) - (cmax * min.x)) / (max.x - min.x));
	r0 = ((rmin - rmax) / (max.y - min.y)) * y0 + (((rmax * max.y) - (rmin * min.y)) / (max.y - min.y));

	// Find the next texel using its vector
	vc0 = patsvec[c0][r0][0];
	vr0 = patsvec[c0][r0][1];

	c1 = c0 + vc0;
	r1 = r0 + vr0;

	// Find the corresponding vertex in the quad space
	x1 = ((c1 * (max.x - min.x)) / (cmax - cmin)) - (((cmin * max.x) - (cmax * min.x)) / (cmax - cmin));
	y1 = ((r1 * (max.y - min.y)) / (rmin - rmax)) - (((cmin * max.x) - (cmax * min.x)) / (cmax - cmin));

	// Use the two vertices to calculate the vector at the vertex
	vx0 = x1 - x0;
	vy0 = y1 - y0;

	icVector3 vxy0(vx0, vy0, vz);

	// Use bilinear interpolation to find the average vector to return
	
	
	//// Current position texels
	//double pc = (cmax * (p.x - min.x)) / (max.x - min.x);
	//double pr = -1 * (rmax * (p.y - max.y)) / (max.y - min.y);

	//vx0 = ((max.x - min.x) / cmax) * patsvec[c0][r0][0] + min.x; // x position of the vector
	//vx1 = ((max.x - min.x) / cmax) * patsvec[c1][r1][0] + min.x; // x position of the vector
	//vx2 = ((max.x - min.x) / cmax) * patsvec[c2][r2][0] + min.x; // x position of the vector
	//vx3 = ((max.x - min.x) / cmax) * patsvec[c3][r3][0] + min.x; // x position of the vector

	//vy0 = -1 * ((max.y - min.y) / rmax) * patsvec[c0][r0][1] + max.y; // y position of the vector
	//vy1 = -1 * ((max.y - min.y) / rmax) * patsvec[c1][r1][1] + max.y; // y position of the vector
	//vy2 = -1 * ((max.y - min.y) / rmax) * patsvec[c2][r2][1] + max.y; // y position of the vector
	//vy3 = -1 * ((max.y - min.y) / rmax) * patsvec[c3][r3][1] + max.y; // y position of the vector

	// The vectors to use in bilinear interpolation
	//icVector3 v11(vx2, vy2, vz);
	//icVector3 v12(vx1, vy1, vz);
	//icVector3 v21(vx3, vy3, vz);
	//icVector3 v22(vx0, vy0, vz);

	//icVector3 v =
	//	(x0 - p.x) / (x0 - x2) * (y0 - p.y) / (y0 - y2) * v11 +
	//	(p.x - x2) / (x0 - x2) * (y0 - p.y) / (y0 - y2) * v21 +
	//	(x0 - p.x) / (x0 - x2) * (p.y - y2) / (y0 - y2) * v12 +
	//	(p.x - x2) / (x0 - x2) * (p.y - y2) / (y0 - y2) * v22;

	//normalize vector
	normalize(vxy0);
	return vxy0;
}


bool quadricRoot(double& r0, double& r1, const double& a, const double& b, const double& c) { //verified
	float m = (b * b) - 4 * a * c;
	if (m < 0) {
		return false;
	} 
	r0 = (-b - std::sqrt(m)) / (2 * a);
	r1 = (-b + std::sqrt(m)) / (2 * a);
	return true;
}

bool singRoot(
	double& r0, double& r1,
	const double& a, const double& b, const double& c, const double& d) { // verified
	double f0 = b - a - (c + d);
	double f1 = (c + d);
	double f2 = a;
	return quadricRoot(r0, r1, f0, f1, f2);
}
