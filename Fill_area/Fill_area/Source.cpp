#include <gl/glut.h>
#include <cmath>
#include <iostream>
using namespace std;

const int MAX_LEVEL = 5;
struct Vec2
{
public:
	float x, y;
	Vec2() {}
	Vec2(float x, float y)
	{
		this->x = x;
		this->y = y;
	}
};

struct Vec2i
{
public:
	int x, y;
	Vec2i(){}
	Vec2i(int x, int y)
	{
		this->x = x;
		this->y = y;
	}
};

void drawDot(Vec2 p)
{
	glBegin(GL_POINTS);
	glVertex2f(p.x, p.y);
	glEnd();
	glFlush();
}

void drawDoti(Vec2i p)
{
	glBegin(GL_POINTS);
	glVertex2i(p.x, p.y);
	glEnd();
	glFlush();
}

void drawLine(Vec2 p1, Vec2 p2)
{
	glBegin(GL_LINES);
	glVertex2f(p1.x, p1.y);
	glVertex2f(p2.x, p2.y);
	glEnd();
	glFlush();
}

Vec2 getMidPoint(Vec2 p1, Vec2 p2, float t = 0.5)
{
	Vec2 p;
	p.x = p1.x * t + p2.x * (1.0 - t);
	p.y = p1.y * t + p2.y * (1.0 - t);
	return p;
}

void draw_bezier(Vec2 p1, Vec2 p2, Vec2 p3, Vec2 p4, int level = 1)
{
	if (level == MAX_LEVEL)
		drawLine(p1, p4);
	else
	{
		Vec2 l1 = p1;
		Vec2 l2 = getMidPoint(p1, p2);
		Vec2 h = getMidPoint(p2, p3);
		Vec2 r3 = getMidPoint(p3, p4);
		Vec2 r4 = p4;
		Vec2 l3 = getMidPoint(l2, h);
		Vec2 r2 = getMidPoint(r3, h);
		Vec2 l4 = getMidPoint(l3, r2);
		Vec2 r1 = l4;
		draw_bezier(l1, l2, l3, l4, level + 1);
		draw_bezier(r1, r2, r3, r4, level + 1);
	}
}

// Fill boundary
GLubyte BoundaryColor[3] = { 255,255,255 };		// ±ßÏßÑÕÉ«£¬ºÚÉ«
GLubyte InteriorColor[3] = { 255,0,0 };			// Ìî³äÑÕÉ«£¬ºìÉ«
GLubyte iPixelColor[3];
int ct = 0;

bool isSameColor(int r1,int g1, int b1, int r2, int g2, int b2)
{
	int dis = 10;	// ÈÝ²î
	if(abs(r1 - r2) <= dis && abs(g1 - g2) <= dis && abs(b1 - b2) <= dis)
		return true;
	else
		return false;
}

void Fill_Boundary_4_Connected(int x, int y)
{
	ct++;
	glReadPixels(x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, &iPixelColor);
	cout << "r: " << (int)iPixelColor[0] << " g: " << (int)iPixelColor[1] << " b: " << (int)iPixelColor[2] << endl;
	if (isSameColor((int)iPixelColor[0], (int)iPixelColor[1], (int)iPixelColor[2],
		(int)BoundaryColor[0], (int)BoundaryColor[1], (int)BoundaryColor[2])
		&& !isSameColor((int)iPixelColor[0], (int)iPixelColor[1], (int)iPixelColor[2],
		(int)InteriorColor[0], (int)InteriorColor[1], (int)InteriorColor[2]))
	{
		drawDoti(Vec2i(x,y));
		Fill_Boundary_4_Connected(x + 1, y);
		Fill_Boundary_4_Connected(x, y - 1);
		Fill_Boundary_4_Connected(x - 1, y);
		Fill_Boundary_4_Connected(x, y + 1);
	}
}

void draw()
{
	glClear(GL_COLOR_BUFFER_BIT);

	glColor3f(0.0, 0.0, 0.0);

	Vec2 p1(80.0, 80.0);
	Vec2 p2(160.0, 400.0);
	Vec2 p3(480.0, 400.0);
	Vec2 p4(560.0, 80.0);
	//	glPointSize(8.0);
	//	drawDot(p1);
	//	drawDot(p2);
	//	drawDot(p3);
	//	drawDot(p4);
	draw_bezier(p1, p2, p3, p4);
	drawLine(p1, p4);
	
	glColor3f(255, 0, 0);
	Fill_Boundary_4_Connected(320, 240);
}

void init()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0, 640.0, 0.0, 480.0);
}

void main(int argc, char ** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(100, 50);
	glutInitWindowSize(640, 480);
	glutCreateWindow("Fill_Boundary");

	init();
	glutDisplayFunc(draw);
	glutMainLoop();
}