/*
 * Instituto Tecnologico de Costa Rica
 * Escuela de Ingenieria en Computacion
 * Computer Graphics
 *
 * Programa: Proyecto 1
 * Archivo:  proyecto_1.c
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>
#include <malloc.h>
#include <GL/glut.h>

#include "proyecto_1.h"

int** findNewCoordinate(int s[][2], int p[][1])
{
    int temp[2][1] = { 0 };

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 1; j++)
            for (int k = 0; k < 2; k++)
                temp[i][j] += (s[i][k] * p[k][j]);

    p[0][0] = temp[0][0];
    p[1][0] = temp[1][0];

    return p;
}

// Scaling the Polygon
//void scale(POLYGON *poly, float sx, float sy)
//{
//    // Initializing the Scaling Matrix.
//    int s[2][2] = {
//            {sx, 0,},
//            {0, sy}
//    };
//
//    int p[2][1];
//
//    // Scaling the triangle
//    for (int i = 0; i < poly->nlines; i++)
//    {
//        p[0][0] = poly[i].lines->p1.x;
//        p[1][0] = poly[i].lines->p1.y;
//        findNewCoordinate(s, p);
//        poly[i].lines->p1.x = p[0][0];
//        poly[i].lines->p1.y = p[1][0];
//
//        p[0][0] = poly[i].lines->p2.x;
//        p[1][0] = poly[i].lines->p2.y;
//        findNewCoordinate(s, p);
//        poly[i].lines->p2.x = p[0][0];
//        poly[i].lines->p2.y = p[1][0];
//    }
//}

//void R(long double alpha) {
//    R_matrix[0][0] = cos(alpha);
//    R_matrix[0][1] = -sin(alpha);
//    R_matrix[0][2] = 0;
//    R_matrix[1][0] = sin(alpha);
//    R_matrix[1][1] = cos(alpha);
//    R_matrix[1][2] = 0;
//    R_matrix[2][0] = 0;
//    R_matrix[2][1] = 0;
//    R_matrix[2][2] = 1;
//}

POINT R(POINT P, int alpha){
    P.x = P.x * cos(alpha) - P.y * sin(alpha);
    P.y = P.x * sin(alpha) + P.y * cos(alpha);
    return P;
}

void R_polygon(POLYGON *poly, int alpha) {
    for (int i = 0; i < poly->nlines; i++) {
        poly->lines[i].p1 = R(poly->lines[i].p1, alpha);
        poly->lines[i].p2 = R(poly->lines[i].p2, alpha);
    }
}

POINT T(POINT P, int Dx, int Dy) {
    P.x += Dx;
    P.y += Dy;
    return P;
}

void T_polygon(POLYGON *poly, long double Dx, long double Dy) {
    for (int i = 0; i < poly->nlines; i++) {
        poly->lines[i].p1 = T(poly->lines[i].p1, Dx, Dy);
        poly->lines[i].p2 = T(poly->lines[i].p2, Dx, Dy);
    }
}

POINT S(POINT P, long double Sx, long double Sy) {
    P.x *= Sx;
    P.y *= Sy;
    return P;
}

void S_polygon(POLYGON *poly, long double Sx, long double Sy) {
    for (int i = 0; i < poly->nlines; i++) {
        poly->lines[i].p1 = S(poly->lines[i].p1, Sx, Sy);
        poly->lines[i].p2 = S(poly->lines[i].p2, Sx, Sy);
    }
}

void plot(int x, int y) {
    buffer[x][y] = global_color;
}

void set_color(double r, double g, double b) {
    global_color.r = r;
    global_color.g = g;
    global_color.b = b;
}

void setLineValues(LINE line) {
    line.b = 0;
    line.m = 0;
    if (line.m != 0) {
        line.delta = -1 / line.m;
    } else {
        line.delta = 0;
    }
}

POLYGON *NewPolygon() {
    POLYGON *np = (POLYGON *) malloc(sizeof(POLYGON));
    np->nlines = -1;
    np->lines = (LINE *) malloc(sizeof(LINE) * 10000); //puntos
}

void AddPolygonLine(POLYGON *poly, int x, int y) {
//    LINE *newLines = (LINE *) realloc(poly->lines, sizeof(LINE)*poly->nlines+1);
//    poly->lines = newLines;
    if (poly->nlines == -1) {
        poly->lines[0].p1.x = x;
        poly->lines[0].p1.y = y;
        poly->nlines++;
    } else if (poly->nlines == 0) {
        poly->lines[0].p2.x = x;
        poly->lines[0].p2.y = y;
        setLineValues(poly->lines[0]);
        poly->nlines++;
    } else {
        poly->lines[poly->nlines].p1.x = poly->lines[poly->nlines - 1].p2.x;
        poly->lines[poly->nlines].p1.y = poly->lines[poly->nlines - 1].p2.y;
        poly->lines[poly->nlines].p2.x = x;
        poly->lines[poly->nlines].p2.y = y;
        setLineValues(poly->lines[poly->nlines]);
        poly->nlines++;
    }
}

void DrawPolygon(POLYGON *poly) {
    if (poly->nlines < 1) { return; }
    for (int i = 0; i < poly->nlines; i++) {
        bresenham(
                poly->lines[i].p1.x,
                poly->lines[i].p1.y,
                poly->lines[i].p2.x,
                poly->lines[i].p2.y
        );
    }
    bresenham(
            poly->lines[poly->nlines - 1].p2.x,
            poly->lines[poly->nlines - 1].p2.y,
            poly->lines[0].p1.x,
            poly->lines[0].p1.y
    );
}

/*
void PaintPolygon(POLYGON * poly){
	if(poly->npoints < 2){return;}
	int initialy = 100000000;
	int finaly = -100000000;
	int highpoint, lowpoint;
	POINT p1,p2;
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

void init() {
    lineCount = 0;
    tool = 1;

    polygons[0] = polyCartago = NewPolygon();
    polygons[1] = polyGuanacaste = NewPolygon();
    polygons[2] = polyHeredia = NewPolygon();
    polygons[3] = polyLimon = NewPolygon();
    polygons[4] = polySanJose = NewPolygon();
    polygons[5] = polyAlajuela = NewPolygon();
    polygons[6] = polyPuntarenas = NewPolygon();
    polygons[7] = polyPuntarenas1 = NewPolygon();
    polygons[8] = poly = NewPolygon();
}



void MyKeyboardFunc(unsigned char Key, int x, int y) {
    switch (Key) {
        case '1':
            tool = 1;
            break;
        case '2':
            tool = 2;
            break;
        case '+':
            S_polygon(poly, 1.1, 1.1);
            break;
        case '-':
            S_polygon(poly, 0.9, 0.9);
            break;
        case 'i':
            R_polygon(poly, 1);
            break;
        case 'j':
            R_polygon(poly, -1);
            break;
    };
    glutPostRedisplay();
}


void SpecialInput(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:
            T_polygon(poly, 0, -10);
            break;
        case GLUT_KEY_DOWN:
            T_polygon(poly, 0, 10);
            break;
        case GLUT_KEY_LEFT:
            T_polygon(poly, -10, 0);
            break;
        case GLUT_KEY_RIGHT:
            T_polygon(poly, 10, 0);
            break;
    };
    glutPostRedisplay();
}


void myMouseFunc(int button, int state, int x, int y) {
    if (tool == 1) {
        if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
            lines[lineCount].p1.x = x;
            lines[lineCount].p1.y = y;
            glutPostRedisplay();
        }

        if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
            lines[lineCount].p2.x = x;
            lines[lineCount].p2.y = y;

            lineCount++;
            glutPostRedisplay();
        }


    }
    if (tool == 2) {
        if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
            AddPolygonLine(poly, x, y);
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

    plot(Xp, Yp);
    while (Xp != x1 || Yp != y1) {
        if (d < 0) {
            Xp += Delta_AX;
            Yp += Delta_AY;
            d += Delta_A;
        } else {
            Xp += Delta_BX;
            Yp += Delta_BY;
            d += Delta_B;
        }
        plot(Xp, Yp);
    }
}

void clear_scene() {
    int i, j;
    for (i = 0; i < H_SIZE; i++) {
        for (j = 0; j < V_SIZE; j++) {
            buffer[i][j].r = 0;
            buffer[i][j].g = 0;
            buffer[i][j].b = 0;
        }
    }
}

void draw_scene() {
    int i, j;

    for (i = 0; i < H_SIZE; i++) {
        for (j = 0; j < V_SIZE; j++) {
            glColor3f(buffer[i][j].r, buffer[i][j].g, buffer[i][j].b);
            glBegin(GL_POINTS);
            glVertex2i(i, j);
            glEnd();
        }
    }
    glFlush();
}

void DrawPolygons(){
    for (int i = 0; i < LEN(polygons); ++i) {
        DrawPolygon(polygons[i]);
    }
}

void renderScene(void) {
    clear_scene();

    for (int i = 0; i < lineCount; ++i) {
        bresenham(lines[i].p1.x, lines[i].p1.y, lines[i].p2.x, lines[i].p2.y);
    }

    DrawPolygons();
    draw_scene();
}

int main(int argc, char **argv) {
    int i, j, length;

    if (argc <= 1) {
        printf("Debes ingresar mas parametros...\n");
        return 1;
    }

    printf("Resolución: %s\n", argv[1]);
    length = strlen(argv[1]);
    for (i = 0; i < length; i++)
        if (!isdigit(argv[1][i])) {
            printf("Ingrese una resulución válida\n");
            return 1;
        }

    // Set Res
    sscanf(argv[1], "%d", &H_SIZE);
    sscanf(argv[1], "%d", &V_SIZE);

    buffer = (COLOR **) malloc(H_SIZE * sizeof(COLOR *));
    for (i = 0; i < H_SIZE; i++) {
        buffer[i] = (COLOR *) malloc(V_SIZE * sizeof(COLOR));
    }

    clear_scene();
    set_color(1, 1, 1);

    init();
    glutInit(&argc, argv);

    // init GLUT and create Window
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(H_SIZE, V_SIZE);
    glutCreateWindow("Proyecto 1");
    glutMouseFunc(myMouseFunc);
    glutKeyboardFunc(MyKeyboardFunc);
    glutSpecialFunc(SpecialInput);

    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluOrtho2D(0.0, H_SIZE, V_SIZE, 0.0);
    // register callbacks
    glutDisplayFunc(renderScene);

//    int m2[2][1] = {1,2};
//    int m1[2][2] = {1,2,3,4};
//
//    int **result = findNewCoordinate(m1, m2);
//    printf("%d ", result[0][0]);

    FILE *file;

    char* rutas[] = {
            "mapas/cartago.txt",
            "mapas/guanacaste.txt",
            "mapas/heredia.txt",
            "mapas/limon.txt",
            "mapas/sanjose.txt",
            "mapas/alajuela.txt",
            "mapas/puntarenas.txt",
            "mapas/puntarenas1.txt",
    };
    for (int k = 0; k <  LEN(rutas); ++k) {
        file = fopen(rutas[k], "r");
        int x;
        int y;

        i = 0;
        char line[4098];
        while (fgets(line, 4098, file))
        {
            // double row[ssParams->nreal + 1];
            char* tmp = strdup(line);

            int j = 0;
            const char* tok;
            for (tok = strtok(line, ","); tok && *tok; j++, tok = strtok(NULL, "\t\n"))
            {
                if(j == 0) x = atoi(tok);
                else y = atoi(tok);
            }


            if(i%1 == 0) {
                printf("%d - %d,%d", k, x, y);
                printf("\n");
                AddPolygonLine(polygons[k], x, y);
            }

            free(tmp);
            i++;
        }
        T_polygon(polygons[k], 50, 75);
        S_polygon(polygons[k], 1.2, 1.2);

    }


    // enter GLUT event processing cycle
    glutMainLoop();

    printf("Fin del programa %s...\n\n", argv[0]);


    return 1;

//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    // register callbacks
//
    return 1;
}
