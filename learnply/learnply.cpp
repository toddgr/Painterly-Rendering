#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <random>

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "trackball.h"
#include "tmatrix.h"

#include "drawUtil.h"

#include "ppm.h"

Polyhedron* poly;
std::vector<PolyLine> lines;
std::vector<icVector3> init_points;
std::vector<icVector3> sources;
std::vector<icVector3> saddles;
std::vector<icVector3> higher_order;
std::vector<icVector3> points;
/*scene related variables*/
const float zoomspeed = 0.9;
int win_width = 1024;
int win_height = 1024;
float aspectRatio = win_width / win_height;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 0.9;

/*
Use keys 1 to 0 to switch among different display modes.
Each display mode can be designed to show one type 
visualization result.

Predefined ones: 
display mode 1: solid rendering
display mode 2: show wireframes
display mode 3: render each quad with colors of vertices
*/
int display_mode = 1;

/*User Interaction related variabes*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = tranlate y, 2 = rotate

// IBFV related variables (Van Wijk 2002)
//https://www.win.tue.nl/~vanwijk/ibfv/
#define NPN		256 //64
#define SCALE	4.0
#define ALPHA	8
float tmax = win_width / (SCALE * NPN);
float dmax = SCALE / win_width;
unsigned char* pixels;

const double STEP = 0.01; // You should experiment to find the optimal step size.
const int STEP_MAX = 10000; // Upper limit of steps to take for tracing each streamline.
std::vector<PolyLine> streamlines; // Used for storing streamlines.

const std::string fname = "../data/image/vader.ppm";
int alpha = (255 * 0.2);
ppm img(fname);
float edge_vectors[NPN][NPN][2]; // For storing the edge field
GLubyte image_colors[NPN][NPN][4];	// For accessing pixel colors
bool blur_image = false;
#define E 2.71828
#define PI 3.1415926

icVector3 min, max;
// min and max texture coords
int rmax = NPN - 1;
int rmin = 0;
int cmin = 0;
int cmax = NPN - 1;

/*****************************************************************************
Global Variables to be messed with for UI
******************************************************************************/
double brush_width = 0.75;
double color_jitter = 0.1;
double brightness = -0.0;
double opacity = 0.5;

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);
void initIBFV();

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void displayIBFV();
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);
void reshape(int width, int height);

/*functions for element picking*/
void display_vertices(GLenum mode, Polyhedron* poly);
void display_quads(GLenum mode, Polyhedron* poly);
void display_selected_vertex(Polyhedron* poly);
void display_selected_quad(Polyhedron* poly);

/*display vis results*/
void display_polyhedron(Polyhedron* poly);
void find_singularities();
void clear_sing_points();
double find_point_x_in_quad(Quad* temp, double x, double y);
double find_point_y_in_quad(Quad* temp, double x, double y);

double sing_prox(icVector2 pos);
Quad* streamline_step(icVector2& cpos, icVector2& npos, Quad* cquad, bool forward);
Vertex* find_vertex(double xx, double yy);
void build_streamline(double x, double y);

void findMinMaxField(icVector3& min, icVector3& max);
icVector3 quadToTexture(double x, double y, double z);
icVector3 findPixelColor(icVector3 v);

// For displaying image
void initImage();
void imageFilter(const std::string& fname);
void displayImage();

// For displaying Sobel filter
void initSobel();
void displaySobel();
void sobelFilter(const std::string& fname);
float gaussFunction(double x, double y, double sigma);
void gaussBlur(const std::string& fname, double sigma);
void initGauss(double std_dev);

void draw_lines(std::vector<icVector3>* points, std::vector<PolyLine>* lines);
void print_test_points();
void print_pixel_color_neighbors(int i, int j);

/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char* argv[])
{
	/*load mesh from ply file*/
	FILE* this_file = fopen("../data/vector_data/v9.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	
	/*initialize the mesh*/
	poly->initialize(); // initialize the mesh
	poly->write_info();


	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Painterly Rendering");

	
	/*initialize openGL*/
	init();

	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	
	/*event processing loop*/
	glutMainLoop();
	
	/*clear memory before exit*/
	poly->finalize();	// finalize everything
	free(pixels);
	clear_sing_points();
	return 0;
}

/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode)
{
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };

	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (aspectRatio >= 1.0) {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, 0.1, 1000);
	}
	else {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, 0.1, 1000);
	}

	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Polyhedron* poly)
{
	glTranslatef(translation[0], translation[1], -3.0);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glEnable(GL_NORMALIZE);
	if (poly->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);

	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* quad = poly->qlist[i];
		quad->singularity = NULL;
	}
}

/******************************************************************************
Initialize IBFV patterns
******************************************************************************/

void initIBFV()
{
	//pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	//memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	//tmax = win_width / (SCALE * NPN);
	//dmax = SCALE / win_width;

	//int lut[256];
	//int phase[NPN][NPN];
	//GLubyte pat[NPN][NPN][4];
	//int i, j, k;

	//for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	//for (i = 0; i < NPN; i++)
	//	for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;

	//for (i = 0; i < NPN; i++)
	//{
	//	for (j = 0; j < NPN; j++)
	//	{
	//		pat[i][j][0] =
	//			pat[i][j][1] =
	//			pat[i][j][2] = lut[(phase[i][j]) % 255];
	//		pat[i][j][3] = ALPHA;
	//	}
	//}
	//
	//glNewList(1, GL_COMPILE);
	//glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
	//glEndList();
}

/******************************************************************************
Pick objects from the scene
******************************************************************************/

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, * ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	ptr = (GLuint*)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	return seed_id;
}

/******************************************************************************
Diaplay all quads for selection
******************************************************************************/

void display_quads(GLenum mode, Polyhedron* this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	for (i = 0; i < this_poly->nquads; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Quad* temp_q = this_poly->qlist[i];
		
		glBegin(GL_POLYGON);
		for (j = 0; j < 4; j++) {
			Vertex* temp_v = temp_q->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

/******************************************************************************
Diaplay all vertices for selection
******************************************************************************/

void display_vertices(GLenum mode, Polyhedron* this_poly)
{
	for (int i = 0; i < this_poly->nverts; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		CHECK_GL_ERROR();

		Vertex* temp_v = this_poly->vlist[i];
		drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15);
	}
	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay selected quad
******************************************************************************/

void display_selected_quad(Polyhedron* this_poly)
{
	if (this_poly->selected_quad == -1)
	{
		return;
	}

	unsigned int i, j;

	glDisable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	Quad* temp_q = this_poly->qlist[this_poly->selected_quad];

	glBegin(GL_POLYGON);
	for (j = 0; j < 4; j++) {
		Vertex* temp_v = temp_q->verts[j];
		glColor3f(1.0, 0.0, 0.0);
		glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
}

/******************************************************************************
Display selected vertex
******************************************************************************/

void display_selected_vertex(Polyhedron* this_poly)
{
	if (this_poly->selected_vertex == -1)
	{
		return;
	}

	Vertex* temp_v = this_poly->vlist[this_poly->selected_vertex];
	drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15, 1.0, 0.0,0.0);

	CHECK_GL_ERROR();
}

/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;
	double sigma = 1.;

	// clear out lines and points
	lines.clear();
	points.clear();

	initImage();	// Initialize image out of input file
	initSobel();
	imageFilter(fname);

	switch (key) {
	case 27:	// set excape key to exit program
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '0':
	{
		display_mode = 7;	// Display mode for original image
		printf("Displaying original image.\n");
		//initImage();	// Initialize image out of input file
		imageFilter(fname);
		glutPostRedisplay();
	}
	break;

	case 'o':
	{
		display_mode = 7;	// Display mode for original image
		printf("Displaying blurred image.\n");
		blur_image = !blur_image;
		//imageFilter(fname);
		initGauss(sigma);
		//gaussBlur(fname, 1.);
		printf("Done.\n");
		glutPostRedisplay();
	}
	break;

	case '1':  // Edge field - Sobel filter implementation
	{
		display_mode = 8;
		printf("Displaying Sobel image.\n");
		imageFilter(fname);
		if (blur_image) {
			initGauss(sigma);
		}
		sobelFilter(fname);
		glutPostRedisplay();
	}
	break;

	case '2':	// streamline display over original image
	{
		display_mode = 3;
		findMinMaxField(min, max);
		std::cout << "\nDrawing streamlines" << std::endl;

		//initSobel();
		if (blur_image) {
			initGauss(sigma);
		}
		
		sobelFilter(fname);

		//for patsvec
		//initImage();
		//imageFilter(fname);

		//if (!streamlines_built) {
			//draw_lines(&points, &streamlines);
			//print_test_points();
			//print_pixel_color_neighbors(64, 64);
			//print_pixel_color_neighbors(140, 140);
			//print_pixel_color_neighbors(192, 192);
			//streamlines_built = true;
		//}

			// make dots along x and y axes
		draw_lines(&points, &streamlines);

		glutPostRedisplay();
	}
	break;

	case '3':	// Brush stroke display
		display_mode = 4;
		findMinMaxField(min, max);
		std::cout << "\nDrawing brush strokes" << std::endl;
		//for patsvec
		//initSobel();
		if (blur_image) {
			initGauss(sigma);
		}
		sobelFilter(fname);

		//if (!streamlines_built) {
			draw_lines(&points, &streamlines);
			//print_test_points();
			//streamlines_built = true;
		//}
		
		glutPostRedisplay();
		break;

	case 's':	// Drawing streamlines and grid dots... test
		display_mode = 4;
		test_lines_drawing(&points, &lines);
		glutPostRedisplay();
		break;

	case 'z':	// solid color display with lighting
		display_mode = 1;
		glutPostRedisplay();
		break;

	case 'x':	// wireframe display
		display_mode = 2;
		glutPostRedisplay();
		break;

	case '5':	// IBFV vector field display
		display_mode = 5;
		saddles.clear();
		higher_order.clear();
		sources.clear();
		find_singularities();
		glutPostRedisplay();
		break;

	case '6':	// add your own display mode
	{
		display_mode = 6;
		saddles.clear();
		higher_order.clear();
		sources.clear();
		find_singularities();
		streamlines.clear();

		for (int k = 0; k < sources.size(); k++)
		{
			icVector3 point = sources[k];
			build_streamline(point.x, point.y);
		}
		for (int k = 0; k < saddles.size(); k++)
		{
			icVector3 point = saddles[k];
			build_streamline(point.x, point.y);
		}
		for (int k = 0; k < higher_order.size(); k++)
		{
			icVector3 point = higher_order[k];
			build_streamline(point.x, point.y);
		}

		glutPostRedisplay();
	}
	break;

	case 'r':	// reset rotation and transformation
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;
	}
}

/******************************************************************************
Callback function for dragging mouse
******************************************************************************/

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case 2:

		Quaternion rvec;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 1:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		
		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0 * x - win_width) / win_width;
			float t = (2.0 * (win_height - y) - win_height) / win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_SHIFT) {  // build up the selection feedback mode

				/*select face*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_quads(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_quad = processHits(hits, selectBuf);
				printf("Selected quad id = %d\n", poly->selected_quad);
				glutPostRedisplay();

				CHECK_GL_ERROR();

			}
			else if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_CTRL)
			{
				/*select vertex*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*  create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_vertices(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_vertex = processHits(hits, selectBuf);
				printf("Selected vert id = %d\n", poly->selected_vertex);

				if (poly->selected_vertex >= 0) {
					Vertex * vtemp = poly->vlist[poly->selected_vertex];
					init_points.push_back(icVector3(vtemp->x, vtemp->y, 0));
				}

				glutPostRedisplay();

			}

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/******************************************************************************
Callback function for window reshaping
******************************************************************************/

void reshape(int width, int height)
{
	win_width = width;
	win_height = height;

	aspectRatio = (float)width / (float)height;

	glViewport(0, 0, width, height);

	set_view(GL_RENDER);

	// reset IBFV pixels buffer
	free(pixels);
	initIBFV();
}

/******************************************************************************
Display IBFV vector field visualization (used for Project 3)
******************************************************************************/

void displayIBFV()
{
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

	glClearColor(0.5, 0.5, 0.5, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw the mesh using pixels and use vector field to advect texture coordinates
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	double modelview_matrix[16], projection_matrix[16];
	int viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glGetIntegerv(GL_VIEWPORT, viewport);

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
			
			icVector2 dp = icVector2(vtemp->vx, vtemp->vy);
			normalize(dp);
			dp *= dmax;

			double dx = -dp.x;
			double dy = -dp.y;

			float px = tx + dx;
			float py = ty + dy;

			glTexCoord2f(px, py);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}

	glEnable(GL_BLEND);

	// blend in noise pattern
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1.0, -1.0, 0.0);
	glScalef(2.0, 2.0, 1.0);

	glCallList(1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// draw the mesh using pixels without advecting texture coords
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
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

/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	set_view(GL_RENDER);
	set_scene(GL_RENDER, poly);

	/*display the mesh*/
	display_polyhedron(poly);

	/*display selected elements*/
	display_selected_vertex(poly);
	display_selected_quad(poly);


	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/

void display_polyhedron(Polyhedron* poly)
{
	unsigned int i, j;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);

	switch (display_mode)
	{
	case 1:	// solid color display with lighting
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.75, 0.75, 0.75, 0.0 };
		GLfloat mat_specular[] = { 0.75, 0.75, 0.75, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 2:	// wireframe display
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(1.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];

			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		glDisable(GL_BLEND);
	}
	break;

	case 3:	// Displays streamlines
	{

		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.75, 0.75, 0.75, 0.0 };
		GLfloat mat_specular[] = { 0.75, 0.75, 0.75, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		//displayImage();

		// draw lines
		for (int k = 0; k < streamlines.size(); k++)
		{
			drawPolyLine(streamlines[k], 1.0, 0.0, 0.0, 0.0);
		}

		glutPostRedisplay();
	}
	break;

	case 4: // brush strokes
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);	
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.0, 0.0, 0.0, 0.0 };
		GLfloat mat_specular[] = { 0.0, 0.0, 0.0, 0.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		initImage();
		imageFilter(fname);
		//displayImage();

		// draw points
		// Here, convert the vertex to pixel space and sample the color at that point on the image.
		// set that color to be that same as the stroke that we are drawing.
		//std::cout << "drawing points now" << std::endl;
		icVector3 prevpoint = icVector3(0., 0., 0.);
		for (int k = 0; k < points.size(); k++)
		{
			icVector3 point = points[k];
			if (point != prevpoint) {
				icVector3 color = findPixelColor(icVector3(point.x, point.y, point.z));
				//std::cout << "drawing point " << k << " at {" << point.x << ", " << point.y << ", " << point.z << "}..." << std::endl;
				drawDot(point.x, point.y, point.z, brush_width, color.x, color.y, color.z, opacity);
			}
			prevpoint = point;
		}
		break;
	}
	break;

	case 5:	// IBFV vector field display
	{
		displayIBFV();
		for (int k = 0; k < sources.size(); k++)
		{
			icVector3 point = sources[k];
			drawDot(point.x, point.y, point.z,  0.15, 1, 0, 0);
		}
		for (int k = 0; k < saddles.size(); k++)
		{
			icVector3 point = saddles[k];
			drawDot(point.x, point.y, point.z, 0.15, 0, 1, 0);
		}
		for (int k = 0; k < higher_order.size(); k++)
		{
			icVector3 point = higher_order[k];
			drawDot(point.x, point.y, point.z, 0.15, 0, 0, 1);
		}
		glutPostRedisplay();
	}
	break;

	case 6: // add your own display mode
	{
		displayIBFV();
		for (int k = 0; k < sources.size(); k++)
		{
			icVector3 point = sources[k];
			drawDot(point.x, point.y, point.z, 0.15, 1, 0, 0);
		}
		for (int k = 0; k < saddles.size(); k++)
		{
			icVector3 point = saddles[k];
			drawDot(point.x, point.y, point.z, 0.15, 0, 1, 0);
		}
		for (int k = 0; k < higher_order.size(); k++)
		{
			icVector3 point = higher_order[k];
			drawDot(point.x, point.y, point.z, 0.15, 0, 0, 1);
		}

		for (int k = 0; k < init_points.size(); k++)
		{
			icVector3 vtemp = init_points[k];
			build_streamline(vtemp.x, vtemp.y);
		}

		for (int k = 0; k < streamlines.size(); k++)
		{
			drawPolyLine(streamlines[k], 1.0, 1.0, 0.0, 0.0);
		}
		glutPostRedisplay();
	}
	break;
	case 7:	// display original image
	{
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
		displayImage();
		glutPostRedisplay();
	}
	break;
	case 8:	// display streamlines
	{
		displaySobel();
		glutPostRedisplay();
	}
	break;

	default:
	{
		// don't draw anything
	}

	}
}


void find_singularities()
{
	// 1. Delete any old singularity points
	clear_sing_points();
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* quad = poly->qlist[i];
		// 2. Find x1, x2, y1, y2, fx1y1, fx2y1, fx1y2, fx2y2, gx1y1, gx2y1, gx1y2, gx2y2
		// x1 is the smallest x coordinate of the 4 vertices
		// x2 is the largest x coordinate of the 4 vertices
		// y1 is the smallest y coordinate of the 4 vertices
		// y2 is the largest y coordinate of the 4 vertices
		// fx1y1 is the vector value x component at the vertex with coordinates (x1,y1)
		// gx1y1 is the vector value y component at the vertex with coordinates (x1,y1)
		// ...
		// to access the vector components of a Vertex* object, use vert->vx and vert->vy
		double x1 = poly->smallest_x(quad);
		double x2 = poly->largest_x(quad);
		double y1 = poly->smallest_y(quad);
		double y2 = poly->largest_y(quad);
		double fx1y1 = find_point_x_in_quad(quad, x1, y1);
		double fx2y1 = find_point_x_in_quad(quad, x2, y1);
		double fx1y2 = find_point_x_in_quad(quad, x1, y2);
		double fx2y2 = find_point_x_in_quad(quad, x2, y2);
		double gx1y1 = find_point_y_in_quad(quad, x1, y1);
		double gx2y1 = find_point_y_in_quad(quad, x2, y1);
		double gx1y2 = find_point_y_in_quad(quad, x1, y2);
		double gx2y2 = find_point_y_in_quad(quad, x2, y2);
		// 3. compute the coefficients for solving the quadratic equation
		double a00 = fx1y1;
		double a10 = fx2y1 - fx1y1;
		double a01 = fx1y2 - fx1y1;
		double a11 = fx1y1 - fx2y1 - fx1y2 + fx2y2;
		double b00 = gx1y1;
		double b10 = gx2y1 - gx1y1;
		double b01 = gx1y2 - gx1y1;
		double b11 = gx1y1 - gx2y1 - gx1y2 + gx2y2;
		double c00 = a11 * b00 - a00 * b11;
		double c10 = a11 * b10 - a10 * b11;
		double c01 = a11 * b01 - a01 * b11;
		// 4. Compute the coefficients of the quadratic equation about s:
		// (-a11*c10)s2 + (-a11*c00 - a01*c10 + a10*c01)s + (a00*c01 - a01*c00) = 0.
		double a = -a11 * c10;
		double b = -a11 * c00 - a01 * c10 + a10 * c01;
		double c = a00 * c01 - a01 * c00;
		// 5. Use the quadratic formula to solve for the s.
		// (check beforehand for complex values or dividing by zero)
		// You will get two values for s because of the ��.
		double s1 = (-b + sqrt(b * b - 4. * a * c)) / 2. / a;
		double s2 = (-b - sqrt(b * b - 4. * a * c)) / 2. / a;
		// 6. Use both values of s to get two corresponding values for t:
		double t1 = -(c00 / c01) - (c10 / c01) * s1;
		double t2 = -(c00 / c01) - (c10 / c01) * s2;
		// 7. For each (s,t) pair, check that both values are between 0 and 1.
		// Either one or none of these pairs will satisfy this condition.
		// If one of the pairs has both components between 0 and 1,
		// then it corresponds to a singularity inside the quad.
		double s, t;

		if ( s1 > 0 && s1 < 1 && t1 > 0 && t1 < 1 ) {
			// (s1, t1)
			s = s1;
			t = t1;
		}
		else if (s2 > 0 && s2 < 1 && t2 > 0 && t2 < 1) {
			// (s2, t2)
			s = s2;
			t = t2;
		}
		else {
			// none of the pairs satisfies
			quad->singularity = NULL;
			continue;
		}
		// 8. Compute the coordinates of the singularity inside the quad using s and t
		// use s to interpolate between x1 and x2 (s tells you how far inbetween
		// x1 and x2 the x coordinate is).
		// use t to interpolate between y1 and y2 (t tells you how far inbetween
		// y1 and y2 the y coordinate is).

		double sing_x, sing_y;
		sing_x = x1 + s * (x2 - x1);
		sing_y = y1 + t * (y2 - y1);
		// 9. Insert the singularity into the quad data structure
		// quad->singularity = new icVector2(x,y);
		// You will need to create a new field in the Quad data structure to store
		// the singularity point. The Quad class definition is located in the 
		// polyhedron.h file.
		quad->singularity = new icVector2(sing_x, sing_y);

		// classify types of singularity
		// 10. calculate the jacobian values dfdx, dfdy, dgdx, dgdy
		double dfdx = (-(y2 - sing_y) * fx1y1 + (y2 - sing_y) * fx2y1 - (sing_y - y1) * fx1y2 + (sing_y - y1) * fx2y2) / (x2 - x1) / (y2 - y1);
		double dfdy = (-(x2 - sing_x) * fx1y1 - (sing_x - x1) * fx2y1 + (x2 - sing_x) * fx1y2 + (sing_x - x1) * fx2y2) / (x2 - x1) / (y2 - y1);
		double dgdx = (-(y2 - sing_y) * gx1y1 + (y2 - sing_y) * gx2y1 - (sing_y - y1) * gx1y2 + (sing_y - y1) * gx2y2) / (x2 - x1) / (y2 - y1);
		double dgdy = (-(x2 - sing_x) * gx1y1 - (sing_x - x1) * gx2y1 + (x2 - sing_x) * gx1y2 + (sing_x - x1) * gx2y2) / (x2 - x1) / (y2 - y1);
		icMatrix2x2 m = icMatrix2x2(dfdx, dfdy, dgdx, dgdy);
		double determ = determinant(m);
		if (determ > 0)
			sources.push_back(icVector3(sing_x, sing_y, 0));
		else if (determ == 0)
			saddles.push_back(icVector3(sing_x, sing_y, 0));
		else
			higher_order.push_back(icVector3(sing_x, sing_y, 0));
	}
}

double find_point_x_in_quad(Quad* temp, double x, double y) {
	for (int i = 0; i < 4; i++) {
		Vertex* vtemp = temp->verts[i];
		if (vtemp->x == x && vtemp->y == y) {
			return vtemp->vx;
		}
	}
}

double find_point_y_in_quad(Quad* temp, double x, double y) {
	for (int i = 0; i < 4; i++) {
		Vertex* vtemp = temp->verts[i];
		if (vtemp->x == x && vtemp->y == y) {
			return vtemp->vy;
		}
	}
}

void clear_sing_points() {
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* quad = poly->qlist[i];
		if (quad->singularity != NULL) {
			// free memory
			delete quad->singularity;
		}
		quad->singularity = NULL;
	}
	return;
}

double sing_prox(icVector2 pos)
{
	double prox = DBL_MAX;
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* quad = poly->qlist[i];
		if (quad->singularity == NULL) continue;
		icVector2 spos = *(quad->singularity);
		// get the (x,y) coordinates of the
		// singularity in an icVector2 object.
		double dist = length(pos - spos);
		if (dist < prox)
			prox = dist;
	}
	return prox;
}

// Find the vertex in the list of vertices in the polyhedron
// This is where we will define the vx and vy of the vertex to be that
// of the edge field
Vertex* find_vertex(double xx, double yy) {
	for (int i = 0; i < poly->nverts; i++) {
		Vertex* v = poly->vlist[i];
		if (std::abs(v->x - xx) < 1.0e-6 && std::abs(v->y - yy) < 1.0e-6) {
			v->vx = quadToTexture(xx, yy, 0).x;
			v->vy = quadToTexture(xx, yy, 0).y;
			return v;
		}
	}
	return NULL;
}

Quad* streamline_step(icVector2& cpos, icVector2& npos, Quad* cquad, bool forward)
{
	//double x1, y1, x2, y2, f11, f12, f21, f22, g11, g21, g12, g22;
	//double f11, f12, f21, f22, g11, g21, g12, g22;
	const double z = 0.; // We don't need anything in the z direction
	Vertex* v11, * v12, * v21, * v22;

	const double x1 = poly->smallest_x(cquad);
	const double x2 = poly->largest_x(cquad);
	const double y1 = poly->smallest_y(cquad);
	const double y2 = poly->largest_y(cquad);

	const icVector3 vxy11 = quadToTexture(x1, y1, 0);
	const icVector3 vxy12 = quadToTexture(x1, y2, 0);
	const icVector3 vxy21 = quadToTexture(x2, y1, 0);
	const icVector3 vxy22 = quadToTexture(x2, y2, 0);

	for (int i = 0; i < poly->nverts; i++) {
		Vertex* v = poly->vlist[i];

		//x1,y1
		if (fabs(v->x - x1) < 1.0e-6 && fabs(v->y - y1) < 1.0e-6) {
			v->vx = vxy11.x;
			v->vy = vxy11.y;
			v11 = v;
			continue;
		}

		//x1,y2
		if (fabs(v->x - x1) < 1.0e-6 && fabs(v->y - y2) < 1.0e-6) {
			v->vx = vxy12.x;
			v->vy = vxy12.y;
			v12 = v;
			continue;
		}

		//x2,y1
		if (fabs(v->x - x2) < 1.0e-6 && fabs(v->y - y1) < 1.0e-6) {
			v->vx = vxy21.x;
			v->vy = vxy21.y;
			v21 = v;
			continue;
		}

		//x2,y2
		if (fabs(v->x - x2) < 1.0e-6 && fabs(v->y - y2) < 1.0e-6) {
			v->vx = vxy22.x;
			v->vy = vxy22.y;
			v22 = v;
			continue;
		}
	}

	//std::cout << "vs found." << std::endl;

	//v11 = find_vertex(x1, y1);
	const double f11 = v11->vx;
	const double g11 = v11->vy;

	//v12 = find_vertex(x1, y2);
	const double f12 = v12->vx;
	const double g12 = v12->vy;

	//v21 = find_vertex(x2, y1);
	const double f21 = v21->vx;
	const double g21 = v21->vy;

	//v22 = find_vertex(x2, y2);
	const double f22 = v22->vx;
	const double g22 = v22->vy;

	const double x0 = cpos.x;
	const double  y0 = cpos.y;
	icVector2 vect;

	const double m1 = (x2 - x0) * (y2 - y0) / (x2 - x1) / (y2 - y1);
	const double m2 = (x0 - x1) * (y2 - y0) / (x2 - x1) / (y2 - y1);
	const double m3 = (x2 - x0) * (y0 - y1) / (x2 - x1) / (y2 - y1);
	const double m4 = (x0 - x1) * (y0 - y1) / (x2 - x1) / (y2 - y1);
	vect.x = m1 * f11 + m2 * f21 + m3 * f12 + m4 * f22;
	vect.y = m1 * g11 + m2 * g21 + m3 * g12 + m4 * g22;
	normalize(vect);
	if (!forward) { vect *= -1.0; }
	npos.x = cpos.x + STEP * vect.x;
	npos.y = cpos.y + STEP * vect.y;

	Quad* nquad = cquad; //guess that the next quad will be the same
	/*check if npos is outside the current quad*/
	if (npos.x < x1 || npos.x > x2 || npos.y < y1 || npos.y > y2)
	{
		// set up local variables
		icVector2 cross_x1, cross_y1, cross_x2, cross_y2;
		double dprod_x1, dprod_y1, dprod_x2, dprod_y2;
		Edge* cross_edge;

		cross_y1 = icVector2(x0 + ((y1 - y0) / vect.y) * vect.x, y1);
		cross_y2 = icVector2(x0 + ((y2 - y0) / vect.y) * vect.x, y2);
		cross_x1 = icVector2(x1, y0 + ((x1 - x0) / vect.x) * vect.y);
		cross_x2 = icVector2(x2, y0 + ((x2 - x0) / vect.x) * vect.y);

		dprod_y1 = dot(vect, cross_y1 - cpos);
		dprod_y2 = dot(vect, cross_y2 - cpos);
		dprod_x1 = dot(vect, cross_x1 - cpos);
		dprod_x2 = dot(vect, cross_x2 - cpos);
		// check y1
		if (cross_y1.x >= x1 && cross_y1.x <= x2 && dprod_y1 > 0 ) {
			npos = cross_y1;
			cross_edge = poly->find_edge(v11, v21);
			nquad = poly->other_quad(cross_edge, cquad);
		}
		//check y2
		else if (cross_y2.x >= x1 && cross_y2.x <= x2 && dprod_y2 > 0) {
			npos = cross_y2;
			cross_edge = poly->find_edge(v12, v22);
			nquad = poly->other_quad(cross_edge, cquad);
		}
		//check x1
		else if (cross_x1.y >= y1 && cross_x1.y <= y2 && dprod_x1 > 0) {
			npos = cross_x1;
			cross_edge = poly->find_edge(v11, v12);
			nquad = poly->other_quad(cross_edge, cquad);
		}
		//check x2
		else if (cross_x2.y >= y1 && cross_x2.y <= y2 && dprod_x2 > 0) {
			npos = cross_x2;
			cross_edge = poly->find_edge(v21, v22);
			nquad = poly->other_quad(cross_edge, cquad);
		}
		// none of the crossing points meets the conditions
		else {
			nquad = poly->find_quad(npos.x, npos.y);
		}

		double proximity = sing_prox(npos);
		if (proximity < STEP) {
			nquad = NULL;
		}
	}
	return nquad;
}

// takes x and y coordinates of a seed point and draws a streamline through that point
void build_streamline(double x, double y)
{
	{
		Quad* cquad = poly->find_quad(x, y);
		icVector2 cpos = icVector2(x, y);
		icVector2 npos;
		int step_counter = 0;
		PolyLine pline;
		while (cquad != NULL && step_counter < STEP_MAX)
		{
			cquad = streamline_step(cpos, npos, cquad, true);
			LineSegment linear_seg = LineSegment(cpos.x, cpos.y, 0, npos.x, npos.y, 0);
			pline.push_back(linear_seg);

			points.push_back(icVector3(cpos.x, cpos.y, 0));	

			cpos = npos;
			step_counter++;
		}
		streamlines.push_back(pline);
	}

	{
		Quad* cquad = poly->find_quad(x, y);
		icVector2 cpos = icVector2(x, y);
		icVector2 npos;
		int step_counter = 0;
		PolyLine pline;
		while (cquad != NULL && step_counter < STEP_MAX)
		{
			cquad = streamline_step(cpos, npos, cquad, false);
			LineSegment linear_seg = LineSegment(cpos.x, cpos.y, 0, npos.x, npos.y, 0);
			pline.push_back(linear_seg);

			// Add a point here, let it sample the color of the image
			points.push_back(icVector3(cpos.x, cpos.y, 0));

			cpos = npos;
			step_counter++;
		}
		streamlines.push_back(pline);
	}
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

// Get c and r values from x and y
icVector3 quadToTexture(double x, double y, double z) {

	// Convert x,y to texture space
	double c = ((y * (cmax - cmin)) +
		((cmin * max.x) - (cmax * min.x))) / (max.x - min.x);
	/*double r = ((y * (rmin - rmax)) +
		((rmax * max.y) - (rmin * min.y))) / (max.y - min.y);*/
	double r = ((x * (rmax - rmin)) +
		((rmin * max.y) - (rmax * min.y))) / (max.y - min.y);
	
	// Find c,r vector in texture space
	int ir = (int)r;
	int ic = (int)c;

	double vc = 0, vr = 0;

	// Make sure that the vector is not out of bounds
	if (ic > 0 || ir > 0 || ic < NPN-1 || ir < NPN-1) {
		vc = edge_vectors[ic][ir][1];
		vr = edge_vectors[ic][ir][0];
	}

	return icVector3(vc, vr, 0.);
}

// Find the color of a pixel given a coordinate on the mesh
icVector3 findPixelColor(icVector3 v) {
	double r = 0., g = 0., b = 0.;
	// v is a vertex in the mesh space
	// convert to pixel space
	int col = ((v.y * (cmax - cmin)) +
		((cmin * max.x) - (cmax * min.x))) / (max.x - min.x);
	int row = ((v.x * (rmax - rmin)) +
		((rmin * max.y) - (rmax * min.y))) / (max.y - min.y);
	
	// find rgb values at a given pixel
	r = image_colors[col][row][0];
	g = image_colors[col][row][1];
	b = image_colors[col][row][2];

	r /= 255;
	g /= 255;
	b /= 255;

	float jitter = -color_jitter + static_cast<float>(rand()) * static_cast<float>(color_jitter + color_jitter) / RAND_MAX;
	/*float red_jitter = -color_jitter + static_cast<float>(rand()) * static_cast<float>(color_jitter + color_jitter) / RAND_MAX;
	float green_jitter = -color_jitter + static_cast<float>(rand()) * static_cast<float>(color_jitter + color_jitter) / RAND_MAX;
	float blue_jitter = -color_jitter + static_cast<float>(rand()) * static_cast<float>(color_jitter + color_jitter) / RAND_MAX;
	*/
	r += jitter + brightness;
	g += jitter + brightness;
	b += jitter + brightness;

	// return
	return icVector3(r, g, b);
}

// Initialize image to be displayed
void initImage()
{
	pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	tmax = win_width / (SCALE * NPN);
	dmax = SCALE / win_width;

	int lut[256];
	int phase[NPN][NPN];
	int i, j, k;

	
	for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	for (i = 0; i < NPN; i++)
		//for (j = 0; j < npn; j++) phase[i][j] = rand() % 256;
		for (j = 0; j < NPN; j++) phase[i][j] = 0.;

	for (i = 0; i < NPN; i++)
	{
		for (j = 0; j < NPN; j++)
		{
			image_colors[i][j][0] =
				image_colors[i][j][1] =
				image_colors[i][j][2] = lut[(phase[i][j]) % 255];
			image_colors[i][j][3] = ALPHA;
		}
	}

	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_colors);
	glEndList();
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

	GLubyte pat[NPN][NPN][4];	// image before filter is applied

	// Set color of each pixel
	int i, j;
	for (i = 0; i < NPN; i++) {		// rows
		for (j = 0; j < NPN; j++) { // columns
			
			image_colors[i][j][0] = img.r[(NPN - i - 1) * NPN + j];
			image_colors[i][j][1] = img.g[(NPN - i - 1) * NPN + j];
			image_colors[i][j][2] = img.b[(NPN - i - 1) * NPN + j];
			image_colors[i][j][3] = alpha;
		}
	}

	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, image_colors);
	glEndList();

}

void displayImage() {

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
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// background for rendering color coding and lighting. Anything that is not the image will be white.
	// Not related to edge display.
	glClearColor(1., 1., 1., 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
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
	glShadeModel(GL_SMOOTH);
}

// Functions for displaying the Sobel filter
void initSobel()
{
	pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	tmax = win_width / (SCALE * NPN);
	dmax = SCALE / win_width;

	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k;

	for (i = 0; i < NPN; i++) lut[i] = i < 127 ? 0 : (NPN - 1);
	for (i = 0; i < NPN; i++)
		//for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;
		for (j = 0; j < NPN; j++) phase[i][j] = 0.;

	for (i = 0; i < NPN; i++)
	{
		for (j = 0; j < NPN; j++)
		{
			pat[i][j][0] =
				pat[i][j][1] =
				pat[i][j][2] = lut[(phase[i][j]) % (NPN - 1)];
			pat[i][j][3] = ALPHA;
		}
	}

	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
	glEndList();
}

void displaySobel() {

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

	//glEnable(GL_BLEND);

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
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
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
	//glDisable(GL_BLEND);
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

	GLubyte intensity[NPN][NPN][4];	// image before filter is applied - intensity
	GLubyte sobelpat[NPN][NPN][4];

	// Set color of each pixel to its intensity
	int i, j;
	for (i = 0; i < NPN; i++) {		// rows
		for (j = 0; j < NPN; j++) { // columns
			//float c = 0.299 * (float)img.r[(NPN - i - 1) * NPN + j]+
			//	0.587 * (float)img.g[(NPN - i - 1) * NPN + j] +
			//	0.114 * (float)img.b[(NPN - i - 1) * NPN + j];

			float c = 0.299 * image_colors[i][j][0] +
				0.587 * image_colors[i][j][1] +
				0.114 * image_colors[i][j][2];

			intensity[i][j][0] = c;
			intensity[i][j][1] = c;
			intensity[i][j][2] = c;
			intensity[i][j][3] = alpha;
		}
	}

	// Sobel filter
	for (i = 1; i < NPN-1; i++) {		// row
		for (j = 1; j < NPN - 1; j++) {	// column
			// Initial magnitude for r,g,b (0,1,2) in the x and y directions
			float magx = 0.0;
			float magy = 0.0;

			// Add magnitude weights from neighbors
			//for (int a = 0; a < 3; a++) {
			//	for (int b = 0; b < 3; b++) {
			//		magx += pat0[i - 1 + a][j - 1 + b][0] * kernelx[a][b];
			//		magy += pat0[i - 1 + a][j - 1 + b][0] * kernely[a][b];
			//	}
			//}

			// Gradients are not producing expected values--
			// Trying the messier version
			// xdir	
			// Top row
			if (i > 0 && j > 0) magx += intensity[i-1][j-1][0] * kernelx[0][0];
			if (j > 0) magx += intensity[i][j-1][0] * kernelx[1][0];
			if (i < NPN-1 && j > 0) magx += intensity[i+1][j-1][0] * kernelx[2][0];

			// Middle row
			if (i > 0) magx += intensity[i - 1][j][0] * kernelx[0][1];
			magx += intensity[i][j][0] * kernelx[1][1];
			if (i < NPN-1) magx += intensity[i+1][j][0] * kernelx[2][1];

			// Bottom row
			if (i >0 && j < NPN-1)magx += intensity[i-1][j+1][0] * kernelx[0][2];
			if (j < NPN-1) magx += intensity[i][j+1][0] * kernelx[1][2];
			if (i < NPN-1 && j < NPN-1) magx += intensity[i+1][j+1][0] * kernelx[2][2];


			//// ydir
			// Top row
			if (i > 0 && j > 0) magy += intensity[i - 1][j - 1][0] * kernely[0][0];
			if (j > 0) magy += intensity[i][j - 1][0] * kernely[1][0];
			if (i < NPN - 1 && j > 0) magy += intensity[i + 1][j - 1][0] * kernely[2][0];
			
			// Middle row
			if (i > 0) magy += intensity[i - 1][j][0] * kernely[0][1];
			magy += intensity[i][j][0] * kernely[1][1];
			if (i < NPN - 1) magy += intensity[i + 1][j][0] * kernely[2][1];

			// Bottom row
			if (i > 0 && j < NPN - 1) magy += intensity[i - 1][j + 1][0] * kernely[0][2];
			if (j < NPN - 1) magy += intensity[i][j + 1][0] * kernely[1][2];
			if (i < NPN - 1 && j < NPN - 1) magy += intensity[i + 1][j + 1][0] * kernely[2][2];

			//mag0x /= NPN;
			//mag0y /= NPN;

			magx = magx;
			magx = magx;

			magy = magy;
			magy = magy;

			// Average values for r, g, b
			float v0 = std::sqrt(magx * magx + magy * magy);	// R gradient
			float v1 = std::sqrt(magx * magx + magy * magy);	// G gradient
			float v2 = std::sqrt(magx * magx + magy * magy);	// B gradient

			// Color value will be either black or white
			if (v0 > 255) v0 = 255;
			if (v0 < 0) v0 = 0;
			if (v1 > 255) v1 = 255;
			if (v1 < 0) v1 = 0;
			if (v2 > 255) v2 = 255;
			if (v2 < 0) v2 = 0;

			// Assign RGB values to current pixel for display
			sobelpat[i][j][0] = v0;
			sobelpat[i][j][1] = v1;
			sobelpat[i][j][2] = v2;
			sobelpat[i][j][3] = alpha;



			// Where the vectors are stored
			auto mag = std::sqrt(magx * magx + magy * magy);
			if (mag != 0) {
				edge_vectors[i][j][0] = magx / mag;
				edge_vectors[i][j][1] = magy / mag;
			}
			else {
				edge_vectors[i][j][0] = 0.;
				edge_vectors[i][j][1] = 0.;
			}

			if (std::isinf(edge_vectors[i][j][0]) || std::isinf(edge_vectors[i][j][1]))
			{
				std::cout << "find infinity" << std::endl;
			}
			

		}
	}
	//for (i = 0; i < NPN; i++) {		// row
	//	for (j = 0; j < NPN; j++) {	// column
	//		if (0 == i || 0 == j || NPN - 1 == i || NPN - 1 == j) {
	//			patsvec[i][j][0] = 0;
	//			patsvec[i][j][1] = 0;
	//		}
	//	}
	//}

	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, sobelpat);
	glEndList();

}

// Calculates the transformation to apply to each pixel in the image
float gaussFunction(double x, double y, double sigma) {
	float exponent = -(x * x + y * y) / (2 * sigma * sigma);
	return (1/(2 * PI * sigma * sigma) * pow(E, exponent));
}

// Apply gaussian blur to image before edge detection
// Aims to create a smoother edge field to work from
void gaussBlur(const std::string& fname, double sigma) {

	// Gaussian kernel-- determined based on the size of sigma
	int kernel_size = ceil(sigma * 3) * 2 + 1;

	double* kernel = new double[kernel_size * kernel_size];
	double sum = 0;
	for (int y = -kernel_size / 2; y <= kernel_size / 2; y++) {
		for (int x = -kernel_size / 2; x <= kernel_size / 2; x++) {
			double value = gaussFunction(x, y, sigma);
			kernel[(y + kernel_size / 2) * kernel_size + x + kernel_size / 2] = value;
			sum += value;
		}
	}

	std::cout << "Gaussian kernel: " << std::endl;
	// Normalize the kernel
	for (int i = 0; i < kernel_size * kernel_size; i++) {
		kernel[i] /= sum;
		std::cout << kernel[i] << ",\t";
		if ((i + 1) % kernel_size == 0 ) {
			std::cout << std::endl;
		}
	}

	

	// gauss texture
	GLubyte image[NPN][NPN][4];	// image before filter is applied - intensity
	//GLubyte gausspat[NPN][NPN][4];

	// For each row and column, check each RGB channel
	int i, j, k; //row, column, color channel
	for (i = 0; i < NPN - 1; i++) {		// row
		for (j = 0; j < NPN - 1; j++) {	// column
			//for (k = 0; k < 3; k++) {
				// If the pixel is a boundary pixel, we keep it the same color
					// (This might be where we can implement edge preservation later as well)

			float magr = 0.0;
			float magg = 0.0;
			float magb = 0.0;


			// Convolve the kernels and the images
			for (int ki = -kernel_size / 2; ki <= kernel_size / 2; ki++) {
				for (int kj = -kernel_size / 2; kj <= kernel_size / 2; kj++) {
					int loc = ((kj + (kernel_size / 2)) * kernel_size) + (ki + (kernel_size / 2));
					magr += image_colors[i - 1 + ki][j - 1 + kj][0] * kernel[loc];
					magg += image_colors[i - 1 + ki][j - 1 + kj][1] * kernel[loc];
					magb += image_colors[i - 1 + ki][j - 1 + kj][2] * kernel[loc];
				}
			}

			// Before we assign this, we need to normalize the magnitudes somehow
			image_colors[i][j][0] = magr;
			image_colors[i][j][1] = magg;
			image_colors[i][j][2] = magb;

			//}
		}
	}

	// Free this memory
	delete[] kernel;

}

// Initialize blurred image to be displayed
void initGauss(double std_dev)
{
	pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	tmax = win_width / (SCALE * NPN);
	dmax = SCALE / win_width;

	for (int i = 0; i < NPN; i++) {		// rows
		for (int j = 0; j < NPN; j++) { // columns

			image_colors[i][j][0] = img.r[(NPN - i - 1) * NPN + j];
			image_colors[i][j][1] = img.g[(NPN - i - 1) * NPN + j];
			image_colors[i][j][2] = img.b[(NPN - i - 1) * NPN + j];
			image_colors[i][j][3] = alpha;
		}
	}

	gaussBlur(fname, std_dev);

	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_colors);
	glEndList();
}

// draws streamlines at consistent points on the image
void draw_lines(std::vector<icVector3>* points, std::vector<PolyLine>* lines)
{
	// make dots along x and y axes
	for (int i = -10; i <= 10; i++)
	{
		for (int j = -10; j <= 10; j++) {
			//std::cout << "building streamline[" << i << "][" << j << "]..." << std::endl;
			build_streamline(i, j);
		}
		if (i == -10) {
			std::cout << "X--------------------" << std::endl;
		}
		if (i == -9) {
			std::cout << "XX-------------------" << std::endl;
		}
		if (i == -8) {
			std::cout << "XXX------------------" << std::endl;
		}
		if (i == -7) {
			std::cout << "XXXX-----------------" << std::endl;
		}
		if (i == -6) {
			std::cout << "XXXXX----------------" << std::endl;
		}
		if (i == -5) {
			std::cout << "XXXXXX---------------" << std::endl;
		}
		if (i == -4) {
			std::cout << "XXXXXXX--------------" << std::endl;
		}
		if (i == -3) {
			std::cout << "XXXXXXXX-------------" << std::endl;
		}
		if (i == -2) {
			std::cout << "XXXXXXXXX------------" << std::endl;
		}
		if (i == -1) {
			std::cout << "XXXXXXXXXX-----------" << std::endl;
		}
		if (i == 0) {
			std::cout << "XXXXXXXXXXX----------" << std::endl;
		}
		if (i == 1) {
			std::cout << "XXXXXXXXXXXX---------" << std::endl;
		}
		if (i == 2) {
			std::cout << "XXXXXXXXXXXXX--------" << std::endl;
		}
		if (i == 3) {
			std::cout << "XXXXXXXXXXXXXX-------" << std::endl;
		}
		if (i == 4) {
			std::cout << "XXXXXXXXXXXXXXX------" << std::endl;
		}
		if (i == 5) {
			std::cout << "XXXXXXXXXXXXXXXX-----" << std::endl;
		}
		if (i == 6) {
			std::cout << "XXXXXXXXXXXXXXXXX----" << std::endl;
		}
		if (i == 7) {
			std::cout << "XXXXXXXXXXXXXXXXXX---" << std::endl;
		}
		if (i == 8) {
			std::cout << "XXXXXXXXXXXXXXXXXXX--" << std::endl;
		}
		if (i == 9) {
			std::cout << "XXXXXXXXXXXXXXXXXXXX-" << std::endl;
		}
		if (i == 10) {
			std::cout << "XXXXXXXXXXXXXXXXXXXXX" << std::endl;
		}
	}
}

void print_test_points() {
	// top left
	std::cout << "top left: (" << quadToTexture(-10, 10, 0).x << ", " << quadToTexture(-10, 10, 0).y << ")" << std::endl;
	//std::cout << "top left: (" << patsvec[cmin][rmax][0] << ", " << patsvec[0][0][1] << ")" << std::endl;
	
	// top middle
	std::cout << "top middle: (" << quadToTexture(0, 10, 0).x << ", " << quadToTexture(0, 10, 0).y << ")" << std::endl;
	//std::cout << "top middle: (" << patsvec[128][0][0] << ", " << patsvec[128][0][1] << ")" << std::endl;
	
	// top right
	std::cout << "top right: (" << quadToTexture(10, 10, 0).x << ", " << quadToTexture(10, 10, 0).y << ")" << std::endl;
	//std::cout << "top right: (" << patsvec[255][0][0] << ", " << patsvec[255][0][1] << ")" << std::endl;
	
	// middle left
	std::cout << "middle left: (" << quadToTexture(-10, 0, 0).x << ", " << quadToTexture(-10, 0, 0).y << ")" << std::endl;
	//std::cout << "middle left: (" << patsvec[0][128][0] << ", " << patsvec[0][128][1] << ")" << std::endl;
	
	// middle
	std::cout << "middle: (" << quadToTexture(0, 0, 0).x << ", " << quadToTexture(0, 0, 0).y << ")" << std::endl;
	//std::cout << "middle: (" << patsvec[128][128][0] << ", " << patsvec[128][128][1] << ")" << std::endl;

	// middle right
	std::cout << "middle right: (" << quadToTexture(10, 0, 0).x << ", " << quadToTexture(10, 0, 0).y << ")" << std::endl;
	//std::cout << "middle right: (" << patsvec[255][128][0] << ", " << patsvec[255][128][1] << ")" << std::endl;

	// bottom left
	std::cout << "bottom left: (" << quadToTexture(-10, -10, 0).x << ", " << quadToTexture(-10, -10, 0).y << ")" << std::endl;
	//std::cout << "bottom left: (" << patsvec[0][255][0] << ", " << patsvec[0][255][1] << ")" << std::endl;

	// bottom middle
	std::cout << "bottom middle: (" << quadToTexture(0, -10, 0).x << ", " << quadToTexture(0, -10, 0).y << ")" << std::endl;
	//std::cout << "bottom middle: (" << patsvec[128][255][0] << ", " << patsvec[128][255][1] << ")" << std::endl;

	// bottom right
	std::cout << "bottom right: (" << quadToTexture(10, -10, 0).x << ", " << quadToTexture(10, -10, 0).y << ")" << std::endl;
	//std::cout << "bottom right: (" << patsvec[255][255][0] << ", " << patsvec[255][255][1] << ")" << std::endl;

	std::cout << "---- MORE TEST POINTS ----" << std::endl;
	for (int i = -10; i <= 10; i++) {
		for (int j = -10; j <= 10; j++) {
			std::cout << "(" << i << ", " << j << "): (" << quadToTexture(i, j, 0).x << ", " << quadToTexture(i, j, 0).y << ")" << std::endl;
		}
	}
}

void print_pixel_color_neighbors(int i, int j) {
	ppm img(fname);

	std::cout << "Pixel neighbors for (" << i << ", " << j << "): " << std::endl;

	std::cout << "{" << (float)(img.r[(NPN - i - 2) * NPN + j - 1]) << ", "
		<< (float)(img.g[(NPN - i - 2) * NPN + j - 1]) << ", "
		<< (float)(img.b[(NPN - i - 2) * NPN + j - 1]) << "}"
		<< " {" << (float)(img.r[(NPN - i - 1) * NPN + j - 1]) << ", "
		<< (float)(img.g[(NPN - i - 1) * NPN + j - 1]) << ", "
		<< (float)(img.b[(NPN - i - 1) * NPN + j - 1]) << "}"
		<< " {" << (float)(img.r[(NPN - i) * NPN + j - 1]) << ", "
		<< (float)(img.g[(NPN - i) * NPN + j - 1]) << ", "
		<< (float)(img.b[(NPN - i) * NPN + j - 1]) << "}" << std::endl;

	std::cout << "{" << (float)(img.r[(NPN - i - 2) * NPN + j]) << ", "
		<< (float)(img.g[(NPN - i - 2) * NPN + j]) << ", "
		<< (float)(img.b[(NPN - i - 2) * NPN + j]) << "}"
		<< " {" << (float)(img.r[(NPN - i - 1) * NPN + j]) << ", "
		<< (float)(img.g[(NPN - i - 1) * NPN + j]) << ", "
		<< (float)(img.b[(NPN - i - 1) * NPN + j]) << "}"
		<< " {" << (float)(img.r[(NPN - i) * NPN + j]) << ", "
		<< (float)(img.g[(NPN - i) * NPN + j]) << ", "
		<< (float)(img.b[(NPN - i) * NPN + j]) << "}" << std::endl;

	std::cout << "{" << (float)(img.r[(NPN - i - 2) * NPN + j + 1]) << ", "
		<< (float)(img.g[(NPN - i - 2) * NPN + j + 1]) << ", "
		<< (float)(img.b[(NPN - i - 2) * NPN + j + 1]) << "}"
		<< " {" << (float)(img.r[(NPN - i - 1) * NPN + j + 1]) << ", "
		<< (float)(img.g[(NPN - i - 1) * NPN + j + 1]) << ", "
		<< (float)(img.b[(NPN - i - 1) * NPN + j + 1]) << "}"
		<< " {" << (float)(img.r[(NPN - i) * NPN + j + 1]) << ", "
		<< (float)(img.g[(NPN - i) * NPN + j + 1]) << ", "
		<< (float)(img.b[(NPN - i) * NPN + j + 1]) << "}" << std::endl;
}