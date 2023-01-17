#include "Project1.h"
#include "polyhedron.h"
#include <iostream>
#include <algorithm>
#include "GL/freeglut.h"
extern Polyhedron* poly;

void findMm(double& M, double& m) {
	m = INFINITY; //minimum	
	M = -m; //maximum
	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		if (vertex->scalar < m) {
			m = vertex->scalar;
		}

		if (vertex->scalar > M) {
			M = vertex->scalar;
		}
	}
}


void HSVtoRGB(icVector3&rgb, const icVector3&hsv) {
	double h = hsv.x;
	double s = hsv.y;
	double v = hsv.z;
	double C = s * v;
	double X = C * (1 - abs(fmod(h / 60.0, 2) - 1));
	double m = v - C;
	double r, g, b;
	if (h >= 0 && h < 60) {
		r = C, g = X, b = 0;
	}
	else if (h >= 60 && h < 120) {
		r = X, g = C, b = 0;
	}
	else if (h >= 120 && h < 180) {
		r = 0, g = C, b = X;
	}
	else if (h >= 180 && h < 240) {
		r = 0, g = X, b = C;
	}
	else if (h >= 240 && h < 300) {
		r = X, g = 0, b = C;
	}
	else {
		r = C, g = 0, b = X;
	}
	rgb.x = (r + m);
	rgb.y = (g + m);
	rgb.z = (b + m);
}


void RGBtoHSV(icVector3& hsv, const icVector3& rgb) {
	double r = rgb.x;
	double g = rgb.y;
	double b = rgb.z;
	// h, s, v = hue, saturation, value
	double cmax = std::max(r, std::max(g, b)); // maximum of r, g, b
	double cmin = std::min(r, std::min(g, b));  // minimum of r, g, b
	double diff = cmax - cmin; // diff of cmax and cmin
	double& h = hsv.x;
	double& s = hsv.y;
	double& v = hsv.z;

	// if cmax and cmin are equal then h = 0
	if (cmax == cmin) {
		h = 0;
	}
	else if (cmax == r) {
		h = fmod(60 * ((g - b) / diff) + 360, 360);
	}
	else if (cmax == g) {
		h = fmod(60 * ((b - r) / diff) + 120, 360);
	}
	else if (cmax == b) {
		h = fmod(60 * ((r - g) / diff) + 240, 360);
	}

	// if cmax equal zero
	if (cmax == 0) {
		s = 0;
	}
	else {
		s = (diff / cmax);
	}

	// compute v
	v = cmax;
}


void project1_1a() {
	double m;
	double M;
	findMm(M, m);
	std::cout << "min: " << m << std::endl;
	std::cout << "max: " << M << std::endl;

	std::cout << "grey scale map" << std::endl;
	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		double s_v = vertex->scalar;
		double gray_ = (s_v - m) / (M - m);
		vertex->R = vertex->G = vertex->B = gray_;
	}
	
	glutPostRedisplay();

}

void project1_1b() {
	double m;
	double M;
	findMm(M, m);
	std::cout << "bi-color map" << std::endl;
	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		double s_v = vertex->scalar;
		icVector3 c1(1.0, 0.0, 0.0);
		icVector3 c2(0.0, 0.0, 1.0); 

		double l = (s_v - m) / (M - m);
		double r = (M - s_v) / (M - m);
		icVector3 c = c1 * l + c2 * r;

		vertex->R = c.x;
		vertex->G = c.y;
		vertex->B = c.z;
	}
	glutPostRedisplay();
}


void project1_1c() {
	double m;
	double M;
	findMm(M, m);
	std::cout << "rainbow map" << std::endl;
	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		double s_v = vertex->scalar;
		icVector3 c1(1.0, 0.0, 0.0);
		icVector3 c2(0.0, 0.0, 1.0);
		icVector3 HSVc1, HSVc2;
		RGBtoHSV(HSVc1, c1);
		RGBtoHSV(HSVc2, c2);
		double l = (s_v - m) / (M - m);
		double r = (M - s_v) / (M - m);
		icVector3 HSVc = HSVc1 * l + HSVc2 * r;
		icVector3 RGBc;
		HSVtoRGB(RGBc, HSVc);

		vertex->R = RGBc.x;
		vertex->G = RGBc.y;
		vertex->B = RGBc.z;
	}
	glutPostRedisplay();
}


void project1_2() {
	double m;
	double M;
	findMm(M, m);
	std::cout << "height map" << std::endl;
	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		double s_v = vertex->scalar;
		double l = (s_v - m) / (M - m);
		vertex->z = 2 * l;

		vertex->R = 1.;
		vertex->G = 0.;
		vertex->B = 1.;
	}
	glutPostRedisplay();
}


void project1_3() {
	double m;
	double M;
	findMm(M, m);
	std::cout << "height map" << std::endl;
	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		double s_v = vertex->scalar;
		double l = (s_v - m) / (M - m);
		vertex->z = 2 * l;

		icVector3 c1(1.0, 0.0, 0.0);
		icVector3 c2(0.0, 0.0, 1.0);

		double r = (M - s_v) / (M - m);
		icVector3 c = c1 * l + c2 * r;

		vertex->R = c.x;
		vertex->G = c.y;
		vertex->B = c.z;
	}
	glutPostRedisplay();
}

void project1_3_exaggerated() {
	double m;
	double M;
	findMm(M, m);
	std::cout << "height map" << std::endl;
	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		double s_v = vertex->scalar;
		double l = (s_v - m) / (M - m);
		vertex->z = 10 * l;

		icVector3 c1(1.0, 0.0, 0.0);
		icVector3 c2(0.0, 0.0, 1.0);

		double r = (M - s_v) / (M - m);
		icVector3 c = c1 * l + c2 * r;

		vertex->R = c.x;
		vertex->G = c.y;
		vertex->B = c.z;
	}
	glutPostRedisplay();
}