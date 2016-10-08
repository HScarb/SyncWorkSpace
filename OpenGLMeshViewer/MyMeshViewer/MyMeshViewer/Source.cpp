/*
	MyMeshViewer
*/
#include <iostream>
#include <string>
#include <OVM/OVMCore/IO/MeshIO.h>
#include <OVM/OVMCore/Mesh/THMesh_ArrayKernelT.h>
#include <gl/glut.h>
#define M_PI 3.14159265358979323846

static float c = M_PI / 180.0f;
static int du = 90, oldmy = -1, oldmx = -1;
static float r = 5.0f, h = 0.0f;

typedef OVM::THMesh_ArrayKernelT<> Mesh;

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH, GL_NICEST);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH, GL_NICEST);

	glLoadIdentity();
//	gluLookAt(r*cos(c*du), h, r*sin(c*du), 0, 0, 0, 0, 1, 0); //从视点看远点,y轴方向(0,1,0)是上方向  
	gluLookAt(r*cos(c*du), 0, r*sin(c*du) + 50, 0, 0, 0, 0, 1, 0); //从视点看远点,y轴方向(0,1,0)是上方向  

	glBegin(GL_TRIANGLE_FAN);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 25.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 25.0f);
	glVertex3f(25.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 25.0f, 0.0f);
	glEnd();

	glPopMatrix();
	glutSwapBuffers();
}

void Mouse(int button, int state, int x, int y)		// mouse click
{
	if (state == GLUT_DOWN)	// 第一次鼠标按下时记录在窗口中的坐标
		oldmx = x, oldmy = y;
}

void onMouseMove(int x, int y)
{
	du += x - oldmx;
	h += 0.03f * (y - oldmy);
//	if (h > 1.0f)
//		h = 1.0f;
//	else if (h <= -1.0f)
//		h = -1.0f;
	oldmx = x;
	oldmy = y;
}

void init()
{
	glEnable(GL_DEPTH_TEST);
}

void reshape(int w, int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(75.0f, (float)w / h, 1.0f, 1000.0f);
	glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char * argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(400, 400);
	glutCreateWindow("OpenGL");
	init();
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutIdleFunc(display);  //设置不断调用显示函数  
	glutMouseFunc(Mouse);
	glutMotionFunc(onMouseMove);
	glutMainLoop();
	return 0;
}