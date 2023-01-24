#include "VectorFieldTopology.h"
#include <iostream>
#define M_PI 3.14159265358979323846

constexpr auto EPSILON = 1.0e-5;
constexpr auto STEP = 0.005;
constexpr auto MIN_K = 0.05;

extern Polyhedron* poly;
extern std::vector<POLYLINE> polylines;
extern std::list<Singularity> singularities;

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

// Streamline tracing  // Verified
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
		for (int e_i = 0; e_i < 4; e_i++) {
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

// Get streamline
void streamlineFB(POLYLINE& line, const icVector3& seed, const double& step, bool forward) {  // Verified
	line.m_vertices.push_back(seed);
	Quad* quad = findQuad(seed);
	icVector3 min, max;
	findMinMaxField(min, max);
	icVector3 currPos = seed;
	double coef = 1.0;
	if (!forward) {
		coef = -1.0;
	}
	while (quad != nullptr) {
		icVector3 currVec = getVector(quad, currPos);
		//v:=0
		if (currVec.length() < EPSILON) {
			break;
		}
		//first euler
		icVector3 nextPos = currPos + step * currVec * coef;
		//on boundary
		if (sinp2Boundary(nextPos, min, max)) {
			line.m_vertices.push_back(nextPos);
			break;
		}
		//get the new quad
		Quad* nextQuad = nullptr;
		//nextQuad = findQuad(nextPos);
		streamlineTrace(nextQuad, quad, currPos, currVec * coef, step, min, max);
		//
		//streamlineTrace(nextQuad, quad, currPos, currVec * coef, step, min, max);
		//update quad, pos, vector
		quad = nextQuad;
		currPos = nextPos;
		line.m_vertices.push_back(currPos);
	}
}

void streamline(POLYLINE& line, const icVector3& seed, const double& step) {  // verified
	streamlineFB(line, seed, step);
	POLYLINE line_back;
	streamlineFB(line_back, seed, step, false);
	line.merge(line_back);
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


// Get the vector from vector field by bilinear interpolation
icVector3 getVector(Quad* q, const icVector3& p) {
	double x1 = q->verts[2]->x;
	double x2 = q->verts[0]->x;
	double y1 = q->verts[2]->y;
	double y2 = q->verts[0]->y;

	icVector3 v11(q->verts[2]->vx, q->verts[2]->vy, q->verts[2]->vz);
	icVector3 v12(q->verts[1]->vx, q->verts[1]->vy, q->verts[1]->vz);
	icVector3 v21(q->verts[3]->vx, q->verts[3]->vy, q->verts[3]->vz);
	icVector3 v22(q->verts[0]->vx, q->verts[0]->vy, q->verts[0]->vz);
	icVector3 v =
		(x2 - p.x) / (x2 - x1) * (y2 - p.y) / (y2 - y1) * v11 +
		(p.x - x1) / (x2 - x1) * (y2 - p.y) / (y2 - y1) * v21 +
		(x2 - p.x) / (x2 - x1) * (p.y - y1) / (y2 - y1) * v12 +
		(p.x - x1) / (x2 - x1) * (p.y - y1) / (y2 - y1) * v22;

	//normalize vector
	//normalize(v);
	return v;
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
void extractSingularity() { // verified
	singularities.clear();
	// Go through faces
	for (int i = 0; i < poly->nquads; i++) {
		// (1)x1y2	1---0 (0)  x2y2
		//			|	|
		// (2)x1y1	2---3 (3) x2y1
		icVector3 vecx1y1 = icVector3(poly->qlist[i]->verts[2]->vx, poly->qlist[i]->verts[2]->vy, poly->qlist[i]->verts[2]->vz);
		icVector3 posx1y1 = icVector3(poly->qlist[i]->verts[2]->x, poly->qlist[i]->verts[2]->y, poly->qlist[i]->verts[2]->z);
		icVector3 vecx2y1 = icVector3(poly->qlist[i]->verts[3]->vx, poly->qlist[i]->verts[3]->vy, poly->qlist[i]->verts[3]->vz);
		icVector3 posx2y1 = icVector3(poly->qlist[i]->verts[3]->x, poly->qlist[i]->verts[3]->y, poly->qlist[i]->verts[3]->z);
		icVector3 vecx2y2 = icVector3(poly->qlist[i]->verts[0]->vx, poly->qlist[i]->verts[0]->vy, poly->qlist[i]->verts[0]->vz);
		icVector3 posx2y2 = icVector3(poly->qlist[i]->verts[0]->x, poly->qlist[i]->verts[0]->y, poly->qlist[i]->verts[0]->z);
		icVector3 vecx1y2 = icVector3(poly->qlist[i]->verts[1]->vx, poly->qlist[i]->verts[1]->vy, poly->qlist[i]->verts[1]->vz);
		icVector3 posx1y2 = icVector3(poly->qlist[i]->verts[1]->x, poly->qlist[i]->verts[1]->y, poly->qlist[i]->verts[1]->z);

		icVector3 pt(0.);
		double f_11 = vecx1y1.x;
		double f_12 = vecx1y2.x;
		double f_21 = vecx2y1.x;
		double f_22 = vecx2y2.x;

		double g_11 = vecx1y1.y;
		double g_12 = vecx1y2.y;
		double g_21 = vecx2y1.y;
		double g_22 = vecx2y2.y;

		double a00 = f_11;
		double a10 = f_21 - f_11;
		double a01 = f_12 - f_11;
		double a11 = f_11 - f_21 - f_12 + f_22;
		double b00 = g_11;
		double b10 = g_21 - g_11;
		double b01 = g_12 - g_11;
		double b11 = g_11 - g_21 - g_12 + g_22;
		double c00 = a11 * b00 - a00 * b11;
		double c10 = a11 * b10 - a10 * b11;
		double c01 = a11 * b01 - a01 * b11;

		if (c01 == 0.0) {
			std::cout << "c01 is 0: " << i << std::endl;
			continue;
		}
		double div = c10 / c01;
		double a = -a11 * div;
		double b = a10 - a01 * div;
		double c = a00 - (a01 + a11) * c00 / c01;
		double s[2];
		bool flag = quadricRoot(s[0], s[1], a, b, c);
		if (!flag)
			continue;
		double t[2];
		t[0] = -c00 / c01 - div * s[0];
		t[1] = -c00 / c01 - div * s[1];

		for (int j = 0; j < 2; j++) {
			if (s[j] > -EPSILON && s[j] < 1 + EPSILON &&
				t[j] > -EPSILON && t[j] < 1 + EPSILON) {
				pt.x = s[j];
				pt.y = t[j];
				Singularity point;
				point.p = pt + posx1y1;
				singularities.push_back(point);
				//check
				double r0 = f_11 + pt.x * (f_21 - f_11) + pt.y * (f_12 - f_11) + pt.x * pt.y * (f_11 - f_21 - f_12 + f_22);
				std::cout << r0 << std::endl;
			}
		}


	}
}

void classifySingularity() {
	for (auto& s : singularities) {

		icVector3 posn = s.p;
		Quad* quad = findQuad(posn);

		int R[4] = { 2, 3, 0, 1 };
		const icVector2 min(quad->verts[R[0]]->x, quad->verts[R[0]]->y);
		const icVector2 max(quad->verts[R[2]]->x, quad->verts[R[2]]->y);
		double len_x = max.x - min.x;
		double len_y = max.y - min.y;

		double dfdx = (-1 / len_x) * (max.y - posn.y) * quad->verts[R[0]]->vx
			+ (1 / len_x) * (max.y - posn.y) * quad->verts[R[1]]->vx
			+ (-1 / len_x) * (posn.y - min.y) * quad->verts[R[3]]->vx
			+ (1 / len_x) * (posn.y - min.y) * quad->verts[R[2]]->vx;

		double dfdy = (-1 / len_y) * (max.x - posn.x) * quad->verts[R[0]]->vx
			+ (posn.x - min.x) * (-1 / len_y) * quad->verts[R[1]]->vx
			+ (max.x - posn.x) * (1 / len_y) * quad->verts[R[3]]->vx
			+ (posn.x - min.x) * (1 / len_y) * quad->verts[R[2]]->vx;

		double dgdx = (-1 / len_x) * (max.y - posn.y) * quad->verts[R[0]]->vy
			+ (1 / len_x) * (max.y - posn.y) * quad->verts[R[1]]->vy
			+ (-1 / len_x) * (posn.y - min.y) * quad->verts[R[3]]->vy
			+ (1 / len_x) * (posn.y - min.y) * quad->verts[R[2]]->vy;

		double dgdy = (max.x - posn.x) * (-1 / len_y) * quad->verts[R[0]]->vy
			+ (posn.x - min.x) * (-1 / len_y) * quad->verts[R[1]]->vy
			+ (max.x - posn.x) * (1 / len_y) * quad->verts[R[3]]->vy
			+ (posn.x - min.x) * (1 / len_y) * quad->verts[R[2]]->vy;

		s.jacobi.entry[0][0] = dfdx;
		s.jacobi.entry[0][1] = dfdy;
		s.jacobi.entry[1][0] = dgdx;
		s.jacobi.entry[1][1] = dgdy;

		double tr = dfdx + dgdy;
		double det = dfdx * dgdy - dfdy * dgdx;

		double delta = tr * tr - 4 * det;

		if (delta >= 0) { //?
			double r1 = 0.5 * (tr + sqrt(delta));
			double r2 = 0.5 * (tr - sqrt(delta));
			if (r1 == 0 && r2 == 0)//unknown
			{
				s.type = -1;
			}
			else if (r1 >= 0 && r2 >= 0)//source
			{
				s.type = 0;
				s.rgb = icVector3(1.0, 0.0, 0.0);
			}
			else if (r1 <= 0 && r2 <= 0) //sink
			{
				s.type = 1;
				s.rgb = icVector3(0.0, 0.0, 1.0);
			}
			if (r1 > 0 && r2 < 0 || r1 < 0 && r2 > 0)//saddle
			{
				s.type = 2;
				s.rgb = icVector3(0.0, 1.0, 0.0);
			}
			else {
				s.type = -1;
			}
		}
		else {
			if (tr == 0)//center
			{
				s.type = 3;
				s.rgb = icVector3(0.0, 1.0, 1.0);
			}
			else {
				s.type = 4;
				s.rgb = icVector3(1.0, 1.0, 0.0);
			}//forus
		}
	}
}


void extractSeparatrix() { // verified, change clear() and MIN_K 
	for (auto& s : singularities) {
		if (s.type == 2)//saddle
		{
			double a = s.jacobi.entry[0][0];
			double b = s.jacobi.entry[0][1];
			double c = s.jacobi.entry[1][0];
			double d = s.jacobi.entry[1][1];

			double Yd = (a + d) / 2;
			double Yr = (c - b) / 2;
			double Ys = std::sqrt((a - d) * (a - d) + (b + c) * (b + c)) / 2;
			double theta = std::atan2(b + c, a - d);
			double phi = std::atan(Yr / Ys);
			double a_cos = std::cos(theta / 2);
			double b_sin = -std::sin(theta / 2);
			double c_sin = std::sin(theta / 2);
			double d_cos = std::cos(theta / 2);
			double a_sin_phi = std::sqrt(std::sin(phi + M_PI / 4));
			double a_cos_phi = std::sqrt(std::cos(phi + M_PI / 4));
			icVector3 maj_v(0.0), min_v(0.0);
			maj_v.x = a_cos * (a_sin_phi + a_cos_phi) + b_sin * (a_sin_phi - a_cos_phi);
			maj_v.y = c_sin * (a_sin_phi + a_cos_phi) + d_cos * (a_sin_phi - a_cos_phi);
			min_v.x = a_cos * (a_sin_phi - a_cos_phi) + b_sin * (a_sin_phi + a_cos_phi);
			min_v.x = c_sin * (a_sin_phi - a_cos_phi) + d_cos * (a_sin_phi + a_cos_phi);
			double k = MIN_K / maj_v.length();

			//outgoing p + kv
			POLYLINE separatrix;
			streamlineFB(separatrix, s.p + k * maj_v, STEP);
			separatrix.m_rgb = icVector3(1.0, 0.0, 0.0);
			polylines.push_back(separatrix);
			//outgoing p-kv;
			separatrix.m_vertices.clear();
			streamlineFB(separatrix, s.p - k * maj_v, STEP);
			separatrix.m_rgb = icVector3(1.0, 0.0, 0.0);
			polylines.push_back(separatrix);
			//incoming p + kw;
			k = MIN_K / min_v.length();
			separatrix.m_vertices.clear();
			streamlineFB(separatrix, s.p + k * min_v, STEP, false);
			separatrix.m_rgb = icVector3(0.0, 0.0, 1.0);
			polylines.push_back(separatrix);
			//incoming p - kw;
			separatrix.m_vertices.clear();
			streamlineFB(separatrix, s.p - k * min_v, STEP, false);
			separatrix.m_rgb = icVector3(0.0, 0.0, 1.0);
			polylines.push_back(separatrix);

		}
	}
}


void display_singularities() {  // verified

	CHECK_GL_ERROR();
	for (auto& sing : singularities) {
		GLUquadric* quadric = gluNewQuadric();
		glPushMatrix();
		glTranslated(sing.p.x, sing.p.y, sing.p.z);
		glColor3f(sing.rgb.x, sing.rgb.y, sing.rgb.z);
		gluSphere(quadric, 0.1, 16, 16);
		glPopMatrix();
		gluDeleteQuadric(quadric);
	}
	glDisable(GL_BLEND);
}