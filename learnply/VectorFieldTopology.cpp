#include "VectorFieldTopology.h"
#include "Sobel.h"
#include <iostream>

#define M_PI 3.14159265358979323846
#define NPN 256 // for sobel-- need to change this later

constexpr auto EPSILON = 1.0e-5;
constexpr auto STEP = 0.005;
constexpr auto MIN_K = 0.05;

// min and max texture coords
int rmin = NPN-1;
int rmax = 0;
int cmin = 0;
int cmax = NPN-1; 

extern Polyhedron* poly;
extern std::vector<POLYLINE> polylines;
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


bool isnp2Boundary(icVector3& currPos, const icVector3& min, const icVector3& max) {
	// Checks if the current position is within the boundary

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

	bool inQuad = false;
	while (!inQuad) {

		// Is the current position outside the field?
		// if yes, the next quad is null
		if (currPos.x < min.x || currPos.x > max.x ||
			currPos.y < min.y || currPos.y > max.y) {
			nextQuad = nullptr;
			return;
		}
		else if (insideQuad(currQuad, currPos)) {
			currPos = currPos + currVec * t;
			nextQuad = findQuad(currPos);
			return;
		}

		double t_ = INFINITY;
		Quad* nextQuad_ = nullptr;
		// For each of the quad edges:
		for (int e_i = 0; e_i < 4; e_i++) {		
			Edge* edge = currQuad->edges[e_i];	
			Vertex* v0 = edge->verts[0];
			Vertex* v1 = edge->verts[1];
			double t_temp;

			// If the size of the edge's x values is significantly small
			if (std::abs(v0->x - v1->x) < EPSILON) {
				// Use x values
				t_temp = (v0->x - currPos.x) / currVec.x;
			}
			else { // Use y values
				t_temp = (v0->y - currPos.y) / currVec.y;
			}
			// If t_temp is positive
			if (t_temp > 0 && t_temp < t_) {
				// Get next quad, the one that isn't the current quad
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
			//t_ = t / 1000;
			//t = t - t_;
			//currPos = currPos + currVec * t_;
			//currQuad = findQuad(currPos);
			
			currPos = currPos + currVec * t;
			nextQuad = findQuad(currPos);
			return;
		}
		else {
			if (t_ >= t) {
				inQuad = true;
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
		//std::cout << "{" << currPos.x << ", " << currPos.y << "," << currPos.z <<
		//	"}: <" << currVec.x << ", " << currVec.y << ", " << currVec.z << ">" << std::endl;

		if (currVec.length() < EPSILON) {				// Control the length of the vector
			break;
		}
		//first euler
		icVector3 nextPos = currPos + step * currVec * coef;
		// Is the next position in the boundary?
		if (isnp2Boundary(nextPos, min, max)) {
			line.m_vertices.push_back(nextPos);	// If yes add to the line
			break;
		}
		//get the new quad
		Quad* nextQuad = nullptr;
		// Trace the streamline into the next quad
		// When does the nextQuad become a nullptr here?
		nextQuad = findQuad(nextPos);
		streamlineTrace(nextQuad, quad, currPos, currVec * coef, step, min, max);
		//update quad, pos, vector
		quad = nextQuad;
		currPos = nextPos;
		line.m_vertices.push_back(currPos);
 	}
}


void streamline(POLYLINE& line, const icVector3& seed, const double& step) {
	// Create streamline forward, backward from seed position
	std::cout << "streamline forward" << std::endl;
	streamlineFB(line, seed, step);				// Create streamline forward
	POLYLINE line_back;
	std::cout << "streamline backward" << std::endl;
	streamlineFB(line_back, seed, step, false);	// Create streamline backward
	line.merge(line_back);						// Merge the two together
}


void drawstreamlines(const double step, const double lineMult) {
	POLYLINE line;

	// Find the minimum and maximum coordinates
	findMinMaxField(min, max);
	for (double i = min.x; i < max.x; i = i + (1/lineMult)) {	// Display streamlines based on the number of lines to display (lineMult)
		line.m_vertices.clear();
		streamline(line, icVector3(i, i, 0), step);		// step = 0.001 typically
		line.m_rgb = icVector3(1.0, 0.0, 0.0);			// Streamlines are white for now
		polylines.push_back(line);						// Add line to polylines
		printf("streamline drawn\n");
	}
}

void drawstreamlinestest(const double step, const double lineMult) {
	// Draws the streamlines for given points: corners, edges, and center. 
	// For our sample image of the ukrainian flag, the only streamline that should
	// draw is the absolute middle.
	POLYLINE line;

	// Find the minimum and maximum coordinates
	findMinMaxField(min, max);
	//for (double i = min.x; i < max.x; i = i + (1/lineMult)) {	// Display streamlines based on the number of lines to display (lineMult)
	//	line.m_vertices.clear();
	//	streamline(line, icVector3(i, i, 0), step);		// step = 0.001 typically
	//	line.m_rgb = icVector3(1.0, 0.0, 0.0);			// Streamlines are white for now
	//	polylines.push_back(line);						// Add line to polylines
	//	printf("streamline drawn\n");
	//}
	//
	// TOP RED
	// MIDDLE GREEN
	// BOTTOM BLUE
	// 
	// top left
	std::cout << "top left" << std::endl;
	line.m_vertices.clear();
	streamline(line, icVector3(-10, 10, 0), step);		// step = 0.001 typically
	line.m_rgb = icVector3(1.0, 0.0, 0.0);			// Streamlines are white for now
	polylines.push_back(line);						// Add line to polylines

	// top middle
	std::cout << "top middle" << std::endl;
	line.m_vertices.clear();
	streamline(line, icVector3(0, 10, 0), step);		// step = 0.001 typically
	line.m_rgb = icVector3(1.0, 0.0, 0.0);			// Streamlines are white for now
	polylines.push_back(line);						// Add line to polylines

	// top right
	std::cout << "top right" << std::endl;
	line.m_vertices.clear();
	streamline(line, icVector3(10, 10, 0), step);		// step = 0.001 typically
	line.m_rgb = icVector3(1.0, 0.0, 0.0);			// Streamlines are white for now
	polylines.push_back(line);						// Add line to polylines

	// middle left
	std::cout << "middle left" << std::endl;
	line.m_vertices.clear();
	streamline(line, icVector3(-10, 0, 0), step);		// step = 0.001 typically
	line.m_rgb = icVector3(0.0, 1.0, 0.0);			// Streamlines are white for now
	polylines.push_back(line);						// Add line to polylines

	// middle
	std::cout << "middle" << std::endl;
	line.m_vertices.clear();
	streamline(line, icVector3(0, 0, 0), step);		// step = 0.001 typically
	line.m_rgb = icVector3(0.0, 1.0, 0.0);			// Streamlines are white for now
	polylines.push_back(line);						// Add line to polylines

	// middle right
	std::cout << "middle right" << std::endl;
	line.m_vertices.clear();
	streamline(line, icVector3(10, 0, 0), step);		// step = 0.001 typically
	line.m_rgb = icVector3(0.0, 1.0, 0.0);			// Streamlines are white for now
	polylines.push_back(line);						// Add line to polylines

	// bottom left
	std::cout << "bottom left" << std::endl;
	line.m_vertices.clear();
	streamline(line, icVector3(-10, -10, 0), step);		// step = 0.001 typically
	line.m_rgb = icVector3(0.0, 0.0, 1.0);			// Streamlines are white for now
	polylines.push_back(line);						// Add line to polylines

	// bottom middle
	std::cout << "bottom middle" << std::endl;
	line.m_vertices.clear();
	streamline(line, icVector3(0, -10, 0), step);		// step = 0.001 typically
	line.m_rgb = icVector3(0.0, 0.0, 1.0);			// Streamlines are white for now
	polylines.push_back(line);						// Add line to polylines

	// bottom right
	std::cout << "bottom right" << std::endl;
	line.m_vertices.clear();
	streamline(line, icVector3(10, -10, 0), step);		// step = 0.001 typically
	line.m_rgb = icVector3(0.0, 0.0, 1.0);			// Streamlines are white for now
	polylines.push_back(line);						// Add line to polylines
}

// Find minimum and maximum coordinates for the polyhedron
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
	std::cout << "texture min: {" << cmin << ", " << rmax << ", 0}" << std::endl;
	std::cout << "texture max: {" << cmax << ", " << rmin << ", 0}" << std::endl;

}

Quad* findQuad(const icVector3& v) {
	// Find the next quad
	for (int i = 0; i < poly->nquads; i++) {
		Quad* qtemp = poly->qlist[i];
		if (insideQuad(qtemp, v))
			return qtemp;
	}
	return nullptr;
}

bool insideQuad(const Quad* q, const icVector3& p) {
	// Check if the point is inside the quad
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

// Get the c and r values at a given vertex x,y
icVector3 quadToTexture(double x, double y, double z) {

	double c = (((cmax - cmin) / (max.x - min.x)) * x) +
		(((cmin * max.x) - (cmax * min.x)) / (max.x - min.x));
	double r = ((rmin - rmax) / (max.y - min.y)) * y + 
		(((max.y * rmax) - (min.y * rmin)) / (max.y - min.y));

	return icVector3(c, r, z);
}

// Find the next texel from c,r vector
icVector3 getNextTexel(double c, double r, double w) {
	int ir = (int)r;
	int ic = (int)c;

	double vc = patsvec[ic][ir][0];
	double vr = patsvec[ic][ir][1];

	double cprime = c + vc;
	double rprime = r + vr;

	return icVector3(cprime, rprime, w);
}

// Get the x and y values from a given texel c, r
icVector3 textureToQuad(double r, double c, double w) {

	double x = (((c * (max.x - min.x)) / (cmax - cmin))) -
		(((cmin * max.x) - (cmax * min.x)) / (cmax - cmin));
	double y = (r * ((min.y - max.y) / (rmax - rmin))) +
		(((rmax * max.y) - (rmin * min.y)) / (rmax - rmin));

	return icVector3(x, y, w);
}

// Get the vector from vector field by *bilinear interpolation*
icVector3 getVector(Quad* q, const icVector3& p) {
	// min.x, max.x,
	// min.y, max.y
	// rmin, rmax
	// cmin, cmax
	double x0, x1, x2, x3,
		x0p, x1p, x2p, x3p,
		y0, y1, y2, y3,
		y0p, y1p, y2p, y3p,
		vc0, vc1, vc2, vc3,
		vr0, vr1, vr2, vr3,
		vx0, vx1, vx2, vx3,
		vy0, vy1, vy2, vy3;

	double vz = 0.; // Need to think about what to have for z vector

	int r0, r1, r2, r3,
		r0p, r1p, r2p, r3p,
		c0, c1, c2, c3,
		c0p, c1p, c2p, c3p;

	// Get the vertices from the quad space
	x0 = q->verts[0]->x;
	y0 = q->verts[0]->y;

	x2 = q->verts[2]->x;
	y2 = q->verts[2]->y;

	x1 = q->verts[1]->x;
	y1 = q->verts[1]->y;

	x3 = q->verts[3]->x;
	y3 = q->verts[3]->y;

	// Get the corresponding texels in the texture space
	
	icVector3 cr0 = quadToTexture(x0, y0, vz);
	icVector3 cr1 = quadToTexture(x1, y1, vz);
	icVector3 cr2 = quadToTexture(x2, y2, vz);
	icVector3 cr3 = quadToTexture(x3, y3, vz);

	// Find the next texel using its vector
	// ... might be the part causing issues. Check vector validity
	// Do we need to be finding the vectors of each vertex in the quad?
	icVector3 cr0p = getNextTexel(cr0.x, cr0.y, vz);
	icVector3 cr1p = getNextTexel(cr1.x, cr1.y, vz);
	icVector3 cr2p = getNextTexel(cr2.x, cr2.y, vz);
	icVector3 cr3p = getNextTexel(cr3.x, cr3.y, vz);


	// Find the corresponding vertex in the quad space
	icVector3 xy0p = textureToQuad(cr0p.x, cr0p.y, vz);
	icVector3 xy1p = textureToQuad(cr1p.x, cr1p.y, vz);
	icVector3 xy2p = textureToQuad(cr2p.x, cr2p.y, vz);
	icVector3 xy3p = textureToQuad(cr3p.x, cr3p.y, vz);

	// Use the two vertices to calculate the vector at the vertex
	icVector3 vxy0((xy0p.x - x0), (xy0p.y - y0), vz);
	icVector3 vxy1((xy1p.x - x1), (xy1p.y - y1), vz);
	icVector3 vxy2((xy2p.x - x2), (xy2p.y - y2), vz);
	icVector3 vxy3((xy3p.x - x3), (xy3p.y - y3), vz);

	x2 = vxy2.x;
	y2 = vxy2.y;

	// Use bilinear interpolation to find the average vector to return

	icVector3 v =
		(x0 - p.x) / (x0 - x2) * (y0 - p.y) / (y0 - y2) * vxy2 +
		(p.x - x2) / (x0 - x2) * (y0 - p.y) / (y0 - y2) * vxy3 +
		(x0 - p.x) / (x0 - x2) * (p.y - y2) / (y0 - y2) * vxy1 +
		(p.x - x2) / (x0 - x2) * (p.y - y2) / (y0 - y2) * vxy0;

	//normalize vector
	//normalize(v);
	return v;
}

