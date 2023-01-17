#include "Project2.h"
#include "Project1.h"
#include "polyhedron.h"
#include <iostream>
#include <algorithm>
#include "GL/freeglut.h"
#include "Polyline.h"
#include <utility>
#include <list>
#include <vector>

extern Polyhedron* poly;
extern std::vector<POLYLINE> polylines;
extern std::vector<POLYLINE> critPolylines;
std::vector<Vertex> critPoints;
std::vector<Vertex> maxCritPoints;
std::vector<Vertex> minCritPoints;
std::vector<Vertex> saddleCritPoints;

int n = 1000;
int i_div = 200;

void project2_2a() {
	/* Evenly divide the interval [min, max] into N
	   sub-intervals (you will experiment what N is 
	   optimal for each given dataset) and then extract
	   the contour for each of the sub-values. Color 
	   the contours with the same color.*/

	double M;
	double m;
	findMm(M, m);

	for (int i = m / n; i <= M / n; i++) {
		std::list<POLYLINE> edgei;
		marchingSquare(edgei, *poly, i * i_div);
		std::vector<POLYLINE> polylinei;
		makePolylineFromEdges(polylinei, edgei);
		for (auto& polyline_ : polylinei) {
			//Color the contours with the same color
			polyline_.m_rgb = icVector3(1, 1, 1);
			polylines.push_back(polyline_);
		}
	}
	glutPostOverlayRedisplay();
}


void project2_2b() {
	/* Now color the contours using different colors with
	   the color scheme from project 1 but render the 
	   underlying surface with a solid color.*/

	double M;
	double m;
	findMm(M, m);

	for (int i = m / n; i <= M / n; i++) {
		std::list<POLYLINE> edgei;
		marchingSquare(edgei, *poly, i * i_div);
		std::vector<POLYLINE> polylinei;
		makePolylineFromEdges(polylinei, edgei);
		for (auto& polyline_ : polylinei) {
			// Color the contours with the color schemes from project 1.
			polyline_.m_rgb = icVector3((i*50)/n, (i * 50) /n, (i * 50) /n);
			polylines.push_back(polyline_);
		}
	}
	glutPostRedisplay();
}


void project2_2c() {
	/* Combine contours with height fields.*/

	project1_3_exaggerated();

	double M;
	double m;
	findMm(M, m);

	for (int i = m / n; i <= M / n; i++) {
		std::list<POLYLINE> edgei;
		marchingSquare(edgei, *poly, i * i_div);
		std::vector<POLYLINE> polylinei;
		makePolylineFromEdges(polylinei, edgei);
		for (auto& polyline_ : polylinei) {
			// Color the contours with solid colors for height fields
			polyline_.m_rgb = icVector3(1., 1., 1.);
			polylines.push_back(polyline_);
		}
	}
	glutPostRedisplay();
}


void project2_2d() {
	/* Combine colors/height/contours.*/

	project1_3_exaggerated();

	double M;
	double m;
	findMm(M, m);

	for (int i = m / n; i <= M / n; i++) {
		std::list<POLYLINE> edgei;
		marchingSquare(edgei, *poly, i * i_div);
		std::vector<POLYLINE> polylinei;
		makePolylineFromEdges(polylinei, edgei);
		for (auto& polyline_ : polylinei) {
			// Color the contours with the color schemes from project 1.
			polyline_.m_rgb = icVector3(std::abs(i * 50 / n), std::abs(i * 50 / n), std::abs(i * 50 / n));
			polylines.push_back(polyline_);
		}
	}
	glutPostRedisplay();
}

void display_critical_points(std::vector<Vertex> &critPoints) {
	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	for (auto& critpoint : critPoints) {
		glPointSize(7);
		glBegin(GL_POINTS);
		glVertex3d(critpoint.x, critpoint.y, critpoint.z);
		glEnd();
	}

	glDisable(GL_BLEND);
	glLineWidth(1);
}

void project2_3a() {
	/* Extract and visualize its critical points */

	calculate_critical_points();

	glutPostRedisplay();

}

void project2_3b() {
	/* Extract and visualize all contours containing at least one critical point.*/
	calculate_critical_points();

	// Takes in vector of polylines, vector of vertices
	// outputs 
	// For each vertex in critPoints
		// 

	glutPostOverlayRedisplay();
}

void display_crit_polyline(std::vector<POLYLINE>& polylines) {

}

void calculate_crit_polylines() {

}

void calculate_critical_points() {
	/*Compare a point to its surrounding neighbors
	if it's less than all of them, add to min vector
	if it's greater than all of them, add to max vector*/

	int n = std::sqrt(poly->nverts);

	for (int i = 0; i < poly->nverts; i++) {
		/*check if the selected vertex is an edge
		  in that case, will have less than 8 neighbors*/

		/* n2  -  n1  -  n0
		   |      |      |
		   n4  -  i   -  n3
		   |      |      |
		   n7  -  n6  -  n5  */

		auto vertex = poly->vlist[i];
		auto scalar = vertex->scalar;
		int n0 = i - n - 1;
		int n1 = i - n;
		int n2 = i - n + 1;
		int n3 = i - 1;
		int n4 = i + 1;
		int n5 = i + n - 1;
		int n6 = i + n;
		int n7 = i + n + 1;

		// corners
		if (i == 0) {
			if (scalar > poly->vlist[1]->scalar &&
				scalar > poly->vlist[n]->scalar &&
				scalar > poly->vlist[n + 1]->scalar) {
				maxCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			} else if (scalar < poly->vlist[1]->scalar &&
				scalar < poly->vlist[n]->scalar &&
				scalar < poly->vlist[n + 1]->scalar) {
				minCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
		}
		else if (i == n - 1) { //top left
			if (scalar > poly->vlist[n3]->scalar &&
				scalar > poly->vlist[n5]->scalar &&
				scalar > poly->vlist[n6]->scalar) {
				maxCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
			else if (scalar < poly->vlist[n3]->scalar &&
				scalar < poly->vlist[n5]->scalar &&
				scalar < poly->vlist[n6]->scalar) {
				minCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
		}
		else if (i == (n - 1) * n) { // bottom right
			if (scalar > poly->vlist[n1]->scalar &&
				scalar > poly->vlist[n2]->scalar &&
				scalar > poly->vlist[n4]->scalar) {
				maxCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
			else if (scalar < poly->vlist[n1]->scalar &&
				scalar < poly->vlist[n2]->scalar &&
				scalar < poly->vlist[n4]->scalar) {
				minCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
		}
		else if (i == poly->nverts - 1) { // bottom left
			if (scalar > poly->vlist[n1]->scalar &&
				scalar > poly->vlist[n0]->scalar &&
				scalar > poly->vlist[n3]->scalar) {
				maxCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
			else if (scalar < poly->vlist[n1]->scalar &&
				scalar < poly->vlist[n0]->scalar &&
				scalar < poly->vlist[n3]->scalar) {
				minCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
		}

		//any of these cases can be saddles
		else if (i < n) { // top
			if (scalar > poly->vlist[n3]->scalar &&
				scalar > poly->vlist[n4]->scalar &&
				scalar > poly->vlist[n5]->scalar &&
				scalar > poly->vlist[n6]->scalar &&
				scalar > poly->vlist[n7]->scalar) {
				maxCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			} else if (scalar < poly->vlist[n3]->scalar &&
				scalar < poly->vlist[n4]->scalar &&
				scalar < poly->vlist[n5]->scalar &&
				scalar < poly->vlist[n6]->scalar &&
				scalar < poly->vlist[n7]->scalar) {
				minCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
		}
		else if (i > ((n - 1) * n) && i < poly->nverts - 1) { // bottom
			if (scalar > poly->vlist[n0]->scalar &&
				scalar > poly->vlist[n1]->scalar &&
				scalar > poly->vlist[n2]->scalar &&
				scalar > poly->vlist[n3]->scalar &&
				scalar > poly->vlist[n4]->scalar) {
				maxCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
			else if (scalar < poly->vlist[n0]->scalar &&
				scalar < poly->vlist[n1]->scalar &&
				scalar < poly->vlist[n2]->scalar &&
				scalar < poly->vlist[n3]->scalar &&
				scalar < poly->vlist[n4]->scalar) {
				minCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
		}
		else if (i % n == 0 && i < (n-1) * n) { // right
			if ((scalar > poly->vlist[n1]->scalar) &&
				(scalar > poly->vlist[n2]->scalar) &&
				(scalar > poly->vlist[n4]->scalar) &&
				(scalar > poly->vlist[n6]->scalar) &&
				(scalar > poly->vlist[n7]->scalar)) {
				maxCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			} 
			else if (scalar < poly->vlist[n1]->scalar &&
				scalar < poly->vlist[n2]->scalar &&
				scalar < poly->vlist[n4]->scalar &&
				scalar < poly->vlist[n6]->scalar &&
				scalar < poly->vlist[n7]->scalar) {
				minCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
		}
		else if ((i + 1) % n == 0) { // left
			if (scalar > poly->vlist[n1]->scalar &&
				scalar > poly->vlist[n2]->scalar &&
				scalar > poly->vlist[n4]->scalar &&
				scalar > poly->vlist[n6]->scalar &&
				scalar > poly->vlist[n7]->scalar) {
				maxCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
			else if (scalar < poly->vlist[n1]->scalar &&
				scalar < poly->vlist[n2]->scalar &&
				scalar < poly->vlist[n4]->scalar &&
				scalar < poly->vlist[n6]->scalar &&
				scalar < poly->vlist[n7]->scalar) {
				minCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
		}
		else {  /*compare point i with 8 neighbors*/
			if (scalar > poly->vlist[n0]->scalar &&
				scalar > poly->vlist[n1]->scalar &&
				scalar > poly->vlist[n2]->scalar &&
				scalar > poly->vlist[n3]->scalar &&
				scalar > poly->vlist[n4]->scalar &&
				scalar > poly->vlist[n5]->scalar &&
				scalar > poly->vlist[n6]->scalar &&
				scalar > poly->vlist[n7]->scalar) {
				maxCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			} else if (scalar < poly->vlist[n0]->scalar &&
				scalar < poly->vlist[n1]->scalar &&
				scalar < poly->vlist[n2]->scalar &&
				scalar < poly->vlist[n3]->scalar &&
				scalar < poly->vlist[n4]->scalar &&
				scalar < poly->vlist[n5]->scalar &&
				scalar < poly->vlist[n6]->scalar &&
				scalar < poly->vlist[n7]->scalar) {
				minCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			} // corner neighbors make a saddle
			else if ((scalar < poly->vlist[n7]->scalar &&
				scalar < poly->vlist[n0]->scalar &&
				scalar > poly->vlist[n2]->scalar &&
				scalar > poly->vlist[n5]->scalar ) ||
				(scalar > poly->vlist[n7]->scalar &&
					scalar > poly->vlist[n0]->scalar &&
					scalar < poly->vlist[n2]->scalar &&
					scalar < poly->vlist[n5]->scalar)) {
				saddleCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			} // top/bottom/left/right neighbors make a saddle
			else if ((scalar < poly->vlist[n1]->scalar &&
				scalar < poly->vlist[n6]->scalar &&
				scalar > poly->vlist[n3]->scalar &&
				scalar > poly->vlist[n4]->scalar) ||
				(scalar > poly->vlist[n1]->scalar &&
					scalar > poly->vlist[n6]->scalar &&
					scalar < poly->vlist[n3]->scalar &&
					scalar < poly->vlist[n4]->scalar)) {
				saddleCritPoints.push_back(*vertex);
				critPoints.push_back(*vertex);
			}
		}
	}
}