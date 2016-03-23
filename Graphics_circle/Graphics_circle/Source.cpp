#include <GL/glut.h>
#include <iostream>

class screenPt
{
private:
	GLint x, y;
public:
	screenPt()
	{
		x = y = 0;
	}
	
	void setCoords(GLint xCoordValue, GLint yCoordValue)
	{
		x = xCoordValue;
		y = yCoordValue;
	}

	GLint getx() const
	{
		return x;
	}

	GLint gety() const
	{
		return y;
	}

	void incrementx()
	{
		x++;
	}

	void decrementy()
	{
		y--;
	}
};

void setPixel(GLint xCoord, GLint yCoord)
{
	glBegin(GL_POINTS);
		glVertex2i(xCoord, yCoord);
	glEnd();
}

void circleMidpoint(GLint xc, GLint yc, GLint radius)
{
	screenPt circPt;

	GLint p = 1 - radius;

	circPt.setCoords(0, radius);

	void circlePlotPoints(GLint, GLint, screenPt);
	// Plot the initial point in each circle quadrant.
	circlePlotPoints(xc, yc, circPt);
	// Calculate next point and plot in each octant.
	while(circPt.getx() < circPt.gety())
	{
		circPt.incrementx();
		if (p < 0)
		{
			p += 2 * circPt.getx() + 1;
		}
		else
		{
			circPt.decrementy();
			p += 2 * (circPt.getx() - circPt.gety()) + 1;
		}
		circlePlotPoints(xc, yc, circPt);
	}
}

void circlePlotPoints(GLint xc, GLint yc, screenPt circPt)
{
	setPixel(xc + circPt.getx(), yc + circPt.gety());
	setPixel(xc - circPt.getx(), yc + circPt.gety());
	setPixel(xc + circPt.getx(), yc - circPt.gety());
	setPixel(xc - circPt.getx(), yc - circPt.gety());
	setPixel(xc + circPt.gety(), yc + circPt.getx());
	setPixel(xc - circPt.gety(), yc + circPt.getx());
	setPixel(xc + circPt.gety(), yc - circPt.getx());
	setPixel(xc - circPt.gety(), yc - circPt.getx());
}

void drawCircle()
{
	glClear(GL_COLOR_BUFFER_BIT);           // Clear display window

	glColor3f(1.0, 0.0, 0.0);               // Set line segment color to red
	circleMidpoint(100, 100, 50);

	glFlush();
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);                  // Initialize GLUT
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);    // Set display mode
	glutInitWindowPosition(50, 100);        // Set top-left display-window position
	glutInitWindowSize(400, 300);           // Set display-window width and height
	glutCreateWindow("Circle");      // Create display window

	// init
	glClearColor(1.0, 1.0, 1.0, 1.0);       // set display-window color to white

	glMatrixMode(GL_PROJECTION);            // set projection parameters
	gluOrtho2D(0.0, 200.0, 0.0, 150.0);

	glutDisplayFunc(drawCircle);           // Send graphics to display window
	glutMainLoop();                         // Display everything and wait
}