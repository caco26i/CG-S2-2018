#include <GL/glut.h>
#include <math.h>
#include <stdlib.h>


void bresenham(int, int, int, int);

struct Point {
	int x;
	int y;
};

struct Line {
	struct Point p1;
	struct Point p2;
	int m;
	int b;
	int delta;
};

void setLineValues(struct Line line){
	line.b=0;
	line.m=0;
	if(line.m!=0){
		line.delta = -1/line.m;
	}
	else{
		line.delta = 0;
	}
}

struct Polygon {
	int nlines;
	struct Line* lines;
};

struct Polygon * NewPolygon(){
	struct Polygon * np = (struct Polygon *) malloc(sizeof(struct Polygon));
	np->nlines = -1;
	np->lines = (struct Point*) malloc(50 * sizeof(struct Line));
} 

void AddPolygonLine(struct Polygon * poly, int x, int y){
	if(poly->nlines == -1){
		poly->lines[0].p1.x = x;
		poly->lines[0].p1.y = y;
		poly->nlines++;
	}else if(poly->nlines == 0){
		poly->lines[0].p2.x = x;
		poly->lines[0].p2.y = y;
		setLineValues(poly->lines[0]);
		poly->nlines++;
	}else{
		poly->lines[poly->nlines].p1.x = poly->lines[poly->nlines-1].p2.x;
		poly->lines[poly->nlines].p1.y = poly->lines[poly->nlines-1].p2.y;
		poly->lines[poly->nlines].p2.x = x;
		poly->lines[poly->nlines].p2.y = y;
		setLineValues(poly->lines[poly->nlines]);
		poly->nlines++;
	}
	
}

void DrawPolygon(struct Polygon * poly){
	if(poly->nlines < 1){return;}
	for(int i = 0; i < poly->nlines; i++){
		bresenham(poly->lines[i].p1.x,poly->lines[i].p1.y,poly->lines[i].p2.x,poly->lines[i].p2.y);
	}
	bresenham(poly->lines[poly->nlines-1].p2.x,poly->lines[poly->nlines-1].p2.y,poly->lines[0].p1.x,poly->lines[0].p1.y);
}
/*
void PaintPolygon(struct Polygon * poly){
	if(poly->npoints < 2){return;}
	int initialy = 100000000;
	int finaly = -100000000;
	int highpoint, lowpoint;
	struct Point p1,p2;
	bool secondPoint = false;

	bool lineState[poly->npoints];
	for(int i = 0; i < poly->npoints; i++){
		lineState[i] = false;
		if(poly->points[i].y < initialy) initialy = poly->points[i].y;
		if(poly->points[i].y > finaly) finaly = poly->points[i].y;
	}
	for(int i = initialy ; i <= finaly; i++){
		for(int j = 0; j< poly->npoints; j++){
			if(poly->points[j].y != poly->points[(j+1)%poly->npoints].y){
				highpoint = (poly->points[j].y > poly->points[(j+1)%poly->npoints].y)?j:(j+1)%poly->npoints;
				lowpoint = (poly->points[j].y <= poly->points[(j+1)%poly->npoints].y)?j:(j+1)%poly->npoints;

				if(poly->points[highpoint].y == i){
					lineState[j] = true;
				}
				if(poly->points[lowpoint].y == i-1){
					lineState[j] = false;
				}

				if(lineState[j] == true){
					if(secondPoint){
						p2.x = 0;
						p2.y = i;
						secondPoint = false;

					}else{
						p1.x = 0;
						p1.y = i;
						secondPoint = true;
					}
				}
			}
			
		}
		
	}
	
	bresenham(poly->points[poly->npoints-1].x,poly->points[poly->npoints-1].y,poly->points[0].x,poly->points[0].y);
}
*/
struct Line lines[1000];
int lineCount;
int tool;
struct Polygon* poly;

void init(){
	lineCount = 0;
	tool = 1;
	poly = NewPolygon();
}

void MyKeyboardFunc(unsigned char Key, int x, int y)
{
	switch(Key)
	{
		case '1':tool = 1; break;
		case '2':tool = 2; break;
	};
}

void myMouseFunc(int button, int state, int x, int y)
{
	if(tool == 1){
		if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
			lines[lineCount].p1.x = x;
			lines[lineCount].p1.y = y;

			glutPostRedisplay();
		}

		if(button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
			lines[lineCount].p2.x = x;
			lines[lineCount].p2.y = y;

			lineCount++;
			glutPostRedisplay();
		}
	}
	if(tool == 2){
		if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
			AddPolygonLine(poly, x,y);
			glutPostRedisplay();
		}

	}
	

}

void bresenham(int x0, int y0, int x1, int y1) {
    int Delta_E, Delta_NE, Delta_N, Delta_SE, Delta_S, Delta_A, Delta_B, Delta_AX, Delta_AY, Delta_BX, Delta_BY, d, Xp, Yp, difx, dify, block;

    if (x0 > x1) {
        int aux = x0;
        x0 = x1;
        x1 = aux;
        aux = y0;
        y0 = y1;
        y1 = aux;
    }

    dify = (y1 - y0);
    difx = (x1 - x0);
    Xp = x0;
    Yp = y0;

    if (y0 <= y1) {
        if (abs(difx) >= abs(dify)) {   //Cuadrante 1
            Delta_E = 2 * dify;
    	    Delta_NE = 2 * dify - 2 * difx;
            Delta_AX = 1;
            Delta_AY = 0;
            Delta_BX = 1;
            Delta_BY = 1;
            Delta_A = Delta_E;
            Delta_B = Delta_NE;
            d = 2 * dify - difx;
        } else {                        //Cuadrante 2
            Delta_NE = 2 * dify - 2 * difx;
            Delta_N = -2 * difx;
            Delta_AX = 1;
            Delta_AY = 1;
            Delta_BX = 0;
            Delta_BY = 1;
            Delta_A = Delta_NE;
            Delta_B = Delta_N;
            d = dify - 2 * difx;
        }
    } else {
        if (abs(difx) >= abs(dify)) {   //Cuadrante 8
            Delta_E = 2 * dify;
            Delta_SE = 2 * difx + 2 * dify;
            Delta_BX = 1;
            Delta_BY = 0;
            Delta_AX = 1;
            Delta_AY = -1;
            Delta_B = Delta_E;
            Delta_A = Delta_SE;
            d = 2 * dify + difx;
        } else {                        //Cuadrante 7
            Delta_SE = 2 * difx + 2 * dify;
            Delta_S = 2 * difx;
            Delta_BX = 1;
            Delta_BY = -1;
            Delta_AX = 0;
            Delta_AY = -1;
            Delta_B = Delta_SE;
            Delta_A = Delta_S;
            d = dify + 2 * difx;
        }
    }

    glBegin(GL_POINTS);
	glVertex2i(Xp,Yp);
    while (Xp != x1 || Yp != y1 ) {
        if (d < 0) {
            Xp += Delta_AX;
            Yp += Delta_AY;
            d += Delta_A;
        } else {
            Xp += Delta_BX;
            Yp += Delta_BY;
            d += Delta_B;
        }

        glVertex2i(Xp,Yp);
    }
    glEnd();
}

int max(int a, int b){
	(a>b)?a:b;
}

void renderScene(void) {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glColor3f(1,1,1);

	for (int i = 0; i < lineCount; ++i)
	{
		bresenham(lines[i].p1.x,lines[i].p1.y,lines[i].p2.x,lines[i].p2.y);
	}

	DrawPolygon(poly);

	glutSwapBuffers();


}

int main(int argc, char **argv)
{
	init();
	// init GLUT and create Window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0,0);
	glutInitWindowSize(500,500);
	glutCreateWindow("Lines");
	glutMouseFunc(myMouseFunc);
	glutKeyboardFunc(MyKeyboardFunc);

	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluOrtho2D( 0.0, 500.0, 500.0,0.0 );
	// register callbacks
	glutDisplayFunc(renderScene);


	// enter GLUT event processing cycle
	glutMainLoop();
	
	return 1;

}
