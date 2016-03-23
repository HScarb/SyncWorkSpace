#include <GL/glut.h>
#include <cmath>
const int X0 = 180;
const int Y0 = 15;
const int XEnd = 10;
const int YEnd = 145;
struct Vec2i
{
	int X, Y;
};

void init()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);       // set display-window color to white

	glMatrixMode(GL_PROJECTION);            // set projection parameters
	gluOrtho2D(0.0, 200.0, 0.0, 150.0);
}

void lineDDA()
{
	Vec2i vecs[10000];
	int count = 0;

	int dx = XEnd - X0, dy = YEnd - Y0, steps, k;
	float xIncrement, yIncrement, x = X0, y = Y0;

	if (fabs(dx) > fabs(dy))
		steps = fabs(dx);		// a pixel = a step
	else
		steps = fabs(dy);
	xIncrement = float(dx) / float(steps);
	yIncrement = float(dy) / float(steps);

	vecs[count].X = round(x);
	vecs[count++].Y = round(y);
	for (k = 0; k < steps;k++)
	{
		x += xIncrement;
		y += yIncrement;
		vecs[count].X = round(x);
		vecs[count++].Y = round(y);
	}

	glClear(GL_COLOR_BUFFER_BIT);           // Clear display window

	glColor3f(1.0, 0.0, 0.0);               // Set line segment color to red
	glBegin(GL_POINTS);
	for (k = 0; k < count;k++)
	{
		glVertex2i(vecs[k].X, vecs[k].Y);
	}
	glEnd();

	glFlush();                              // Process all OpenGL routines as quickly as possible
}

// Bresenham line-drawing procedure for |m| < 1.0
void lineBres()
{
	glClear(GL_COLOR_BUFFER_BIT);           // Clear display window

	glColor3f(1.0, 0.0, 0.0);               // Set line segment color to red

	int xEnd = XEnd, yEnd = YEnd, x0 = X0, y0 = Y0;
	int dx = fabs(xEnd - x0), dy = fabs(yEnd - y0);
	int p = 2 * dy - dx;
	int twoDy = 2 * dy, twoDyMinusDx = 2 * (dy - dx);
	int x, y;

	// Determine which endpoint to use as start position.
	if (x0 > xEnd)
	{
		x = xEnd;
		y = yEnd;
		xEnd = x0;
	}
	else
	{
		x = x0;
		y = y0;
	}

	glBegin(GL_POINTS);
	glVertex2i(x, y);

	while (x < xEnd)
	{
		x++;
		if (p < 0)
			p += twoDy;
		else
		{
			y++;
			p += twoDyMinusDx;
		}
		glVertex2i(x, y);
	}

	glEnd();
	glFlush();
}

void lineSegment()
{
	glClear(GL_COLOR_BUFFER_BIT);           // Clear display window

	glColor3f(1.0, 0.0, 0.0);               // Set line segment color to red
	glBegin(GL_LINES);
	glVertex2i(180, 15);
	glVertex2i(10, 145);
	glEnd();

	glFlush();                              // Process all OpenGL routines as quickly as possible
}

void main(int argc, char ** argv)
{
	glutInit(&argc, argv);                  // Initialize GLUT
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);    // Set display mode
	glutInitWindowPosition(50, 100);        // Set top-left display-window position
	glutInitWindowSize(400, 300);           // Set display-window width and height
	glutCreateWindow("An Example OpenGL Program");      // Create display window

	init();                                 // Execute initialization procedure
	glutDisplayFunc(lineBres);           // Send graphics to display window
	glutMainLoop();                         // Display everything and wait
}