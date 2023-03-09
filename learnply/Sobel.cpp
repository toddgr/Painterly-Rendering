#include "Sobel.h"
#include "gl/glew.h"
#include "polyhedron.h"
#include "VectorFieldTopology.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "ppm.h"

extern Polyhedron* poly;
#define NPN	256 //64
#define SCALE	4.0
#define ALPHA	8

extern bool original_image;
extern bool flow_image;
extern int win_width;
extern int win_height;
extern unsigned char* pixels;
extern unsigned char* original_pixels;
extern std::vector<POLYLINE> polylines;

float tmax1 = win_width; // / (SCALE * NPN);
float dmax1 = SCALE / win_width;

GLubyte patsvec[NPN][NPN][2];

int alpha1 = (255 * 0.2);

void initSobel()
{
	pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	tmax1 = win_width / (SCALE * NPN);
	dmax1 = SCALE / win_width;

	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k;

	for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	for (i = 0; i < NPN; i++)
		//for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;
		for (j = 0; j < NPN; j++) phase[i][j] = 0.;

	for (i = 0; i < NPN; i++)
	{
		for (j = 0; j < NPN; j++)
		{
			pat[i][j][0] =
				pat[i][j][1] =
				pat[i][j][2] = lut[(phase[i][j]) % 255];
			pat[i][j][3] = ALPHA;
		}
	}

	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
	glEndList();
}

void displaySobel() {
	// Can we display *just* the Sobel filter without needing to toggle IBFV?
	// Perhaps that is the complication here.

	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_BLEND);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// background for rendering color coding and lighting. If not an edge, pixel will be black.
	glClearColor(0., 0., 0., 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw the mesh using pixels and use vector field to advect texture coordinates
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	double modelview_matrix[16], projection_matrix[16];
	int viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glGetIntegerv(GL_VIEWPORT, viewport);

	glEnable(GL_BLEND);

	// blend in noise pattern 
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1., -1., 0.);
	glScalef(2.0, 2.0, 1.0);

	glCallList(1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax1); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax1, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax1, tmax1); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// background for rendering color coding and lighting. Anything that is not the image will be white.
	// Not related to edge display.
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	// for each quad on the polygon
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* qtemp = poly->qlist[i];
		glBegin(GL_QUADS);
		// for each vertex in a polygon:
		for (int j = 0; j < 4; j++)
		{
			Vertex* vtemp = qtemp->verts[j];
			double tx, ty, dummy;
			gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
				modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);
			tx = tx / win_width;
			ty = ty / win_height;
			glTexCoord2f(tx, ty);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
}


void sobelFilter(const std::string& fname) {

	float kernelx[3][3] = {
		{-1, 0, 1},
		{-2, 0, 2},
		{-1, 0, 1} };

	float kernely[3][3] = {
		{-1, -2, -1},
		{0, 0, 0},
		{1, 2, 1} };

	ppm img(fname);
	GLubyte pat0[NPN][NPN][4];	// image before filter is applied - intensity
	GLubyte pat[NPN][NPN][4];	// image after filter is applied - edge field?
								// add a fifth and sixth element to store the vx and vy elements?

	// Set color of each pixel to its intensity
	int i, j;
	for (i = 0; i < NPN; i++) {		// rows
		for (j = 0; j < NPN; j++) { // columns
			float c = 0.299 * img.r[(NPN - i - 1) * NPN + j] +
				0.587 * img.g[(NPN - i - 1) * NPN + j] +
				0.114 * img.b[(NPN - i - 1) * NPN + j];
			pat0[i][j][0] = c;
			pat0[i][j][1] = c;
			pat0[i][j][2] = c;
			pat0[i][j][3] = alpha1;
		}
	}

	// Sobel filter
	for (i = 1; i < NPN - 1; i++) {		// row
		for (j = 1; j < NPN - 1; j++) {	// column
			// Initial magnitude for r,g,b (0,1,2) in the x and y directions
			float mag0x = 0.0;
			float mag1x = 0.0;
			float mag2x = 0.0;
			float mag0y = 0.0;
			float mag1y = 0.0;
			float mag2y = 0.0;

			for (int a = 0; a < 3; a++) {
				for (int b = 0; b < 3; b++) {
					mag0x += pat0[i - 1 + a][j - 1 + b][0] * kernelx[a][b];
					mag1x += pat0[i - 1 + a][j - 1 + b][1] * kernelx[a][b];
					mag2x += pat0[i - 1 + a][j - 1 + b][2] * kernelx[a][b];

					mag0y += pat0[i - 1 + a][j - 1 + b][0] * kernely[a][b];
					mag1y += pat0[i - 1 + a][j - 1 + b][1] * kernely[a][b];
					mag2y += pat0[i - 1 + a][j - 1 + b][2] * kernely[a][b];
				}
			}

			// Average values for r, g, b
			float v0 = std::sqrt(mag0x * mag0x + mag0y * mag0y);	// R gradient
			float v1 = std::sqrt(mag1x * mag1x + mag1y * mag1y);	// G gradient
			float v2 = std::sqrt(mag2x * mag2x + mag2y * mag2y);	// B gradient

			// Color value will be either black or white
			if (v0 > 255) v0 = 255;
			if (v0 < 0) v0 = 0;
			if (v1 > 255) v1 = 255;
			if (v1 < 0) v1 = 0;
			if (v2 > 255) v2 = 255;
			if (v2 < 0) v2 = 0;

			// Assign RGB values to current pixel for display
			pat[i][j][0] = v0;
			pat[i][j][1] = v1;
			pat[i][j][2] = v2;
			pat[i][j][3] = alpha1;

			//Assign magnitude values for x and y directions:
			patsvec[i][j][0] = mag0x;
			patsvec[i][j][1] = mag0y;
			

			//printf("patsvec[%d][%d][x]: %f \
			//	patsvec[% d][% d][y]: %f\n\n", i, j, patsvec[i][j][0], i, j, patsvec[i][j][1]);
		}
	}

	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, pat);
	glEndList();

}

void createEdgeFieldFromSobel() { // should not be void forever. Take in image and poly file?
	// Returns an edge field
	GLubyte edgefield[NPN][NPN][2];	// store vx, vy value per... texel? pixel? What type should this have

	// do all of the things that the original sobel filter does, 
	// but calculate and store each of the vx and vy values for each 9opixel? texel?) that are used in displaysobel

	// but instead of going to display the sobel filter, we assign the values to the quad directly
	// no texture nonsense



	
	// return edge field
}

void initImage()
{
	pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	tmax1 = win_width / (SCALE * NPN);
	dmax1 = SCALE / win_width;

	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k;

	
	for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	for (i = 0; i < NPN; i++)
		//for (j = 0; j < npn; j++) phase[i][j] = rand() % 256;
		for (j = 0; j < NPN; j++) phase[i][j] = 0.;

	for (i = 0; i < NPN; i++)
	{
		for (j = 0; j < NPN; j++)
		{
			pat[i][j][0] =
				pat[i][j][1] =
				pat[i][j][2] = lut[(phase[i][j]) % 255];
			pat[i][j][3] = ALPHA;
		}
	}

	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
	glEndList();
}

void displayImage() {
	// Can we display *just* the Sobel filter without needing to toggle IBFV?
	// Perhaps that is the complication here.

	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_BLEND);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// background for rendering color coding and lighting. If not an edge, pixel will be blue (for testing purposes).
	glClearColor(0.5, 0.5, 0.5, 1.);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw the mesh using pixels and use vector field to advect texture coordinates
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	double modelview_matrix[16], projection_matrix[16];
	int viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glGetIntegerv(GL_VIEWPORT, viewport);

	glEnable(GL_BLEND);

	// blend in pattern 
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1., -1., 0.);
	glScalef(2.0, 2.0, 1.0);

	glCallList(1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax1); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax1, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax1, tmax1); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// background for rendering color coding and lighting. Anything that is not the image will be white.
	// Not related to edge display.
	//glClearColor(1., 1., 1., 1.0);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* qtemp = poly->qlist[i];
		glBegin(GL_QUADS);
		for (int j = 0; j < 4; j++)
		{
			Vertex* vtemp = qtemp->verts[j];
			double tx, ty, dummy;
			gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
				modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);
			tx = tx / win_width;
			ty = ty / win_height;
			glTexCoord2f(tx, ty);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
}


void imageFilter(const std::string& fname) {

	float kernelx[3][3] = {
		{-1, 0, 1},
		{-2, 0, 2},
		{-1, 0, 1} };

	float kernely[3][3] = {
		{-1, -2, -1},
		{0, 0, 0},
		{1, 2, 1} };

	ppm img(fname);
	GLubyte pat[NPN][NPN][4];	// image before filter is applied - intensity
	GLubyte pat0[NPN][NPN][4];	// image after filter is applied - edge field?

	// Set color of each pixel -- something is wrong here
	int i, j;
	for (i = 0; i < NPN; i++) {		// rows
		for (j = 0; j < NPN; j++) { // columns
			
			pat0[i][j][0] = img.g[(NPN - i - 1) * NPN + j];
			pat0[i][j][1] = img.b[(NPN - i - 1) * NPN + j];
			pat0[i][j][2] = img.r[(NPN - i - 1) * NPN + j];
			pat0[i][j][3] = alpha1;
		}
	}

	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, pat0);
	glEndList();

}