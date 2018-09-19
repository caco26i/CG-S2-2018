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

int Xc = 500;
int Yc = 500;
int r = 1;

float Tx, Ty, alpha, Sx, Sy;

float **GT;
float **PT;
float **LT;
float **AT;
float **CT;
float **SJT;
float **HT;
float **SEA;
int textureId;
int drawType;

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

int compare(const void *a, const void *b) {

    POINT *orderA = (POINT *) a;
    POINT *orderB = (POINT *) b;

    return (orderA->x - orderB->x);
}


POINT R(POINT P, float alpha) {
    int px = P.x;
    int py = P.y;
    P.x = (float) px * (float) cos(alpha) - (float) py * (float) sin(alpha);
    P.y = (float) px * (float) sin(alpha) + (float) py * (float) cos(alpha);
    return P;
}

void R_polygon(POLYGON *poly, float alpha) {
    T_polygon(poly, -Xc, -Yc);

    for (int i = 0; i < poly->nlines; i++) {
        poly->lines[i].p1 = R(poly->lines[i].p1, alpha);
        poly->lines[i].p2 = R(poly->lines[i].p2, alpha);
    }
    T_polygon(poly, Xc, Yc);
}


POINT S(POINT P, float Sx, float Sy) {
    P.x *= Sx;
    P.y *= Sy;
    return P;
}

void S_polygon(POLYGON *poly, float Sx, float Sy) {
    T_polygon(poly, -Xc, -Yc);
    for (int i = 0; i < poly->nlines; i++) {
        poly->lines[i].p1 = S(poly->lines[i].p1, Sx, Sy);
        poly->lines[i].p2 = S(poly->lines[i].p2, Sx, Sy);
    }
    T_polygon(poly, Xc, Yc);
}

void R_db(float alpha) {
    for (int i = 0; i < LEN(polygons);
    ++i) {
        R_polygon(polygons[i], alpha);
    }
}

void T_db(long double Dx, long double Dy) {

    for (int i = 0; i < LEN(polygons);
    ++i) {
        T_polygon(polygons[i], Dx, Dy);
    }
}

void S_db(long double Sx, long double Sy) {
    //printf("Sx %d Sy %d\n",Sx,Sy );
    for (int i = 0; i < LEN(polygons);
    ++i) {
        S_polygon(polygons[i], Sx, Sy);
    }
}

// Returns x-value of point of intersectipn of two
// lines
int x_intersect(int x1, int y1, int x2, int y2,
                int x3, int y3, int x4, int y4) {
    int num = (x1 * y2 - y1 * x2) * (x3 - x4) -
              (x1 - x2) * (x3 * y4 - y3 * x4);
    int den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    return num / den;
}

// Returns y-value of point of intersectipn of
// two lines
int y_intersect(int x1, int y1, int x2, int y2,
                int x3, int y3, int x4, int y4) {
    int num = (x1 * y2 - y1 * x2) * (y3 - y4) -
              (y1 - y2) * (x3 * y4 - y3 * x4);
    int den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    return num / den;
}

// This functions clips all the edges w.r.t one clip
// edge of clipping area
void clip(POLYGON *poly,
          int x1, int y1, int x2, int y2) {
    int MAX_POINTS = 1;//HARDCODE
    int new_points[poly->nlines][2];
    int new_poly_size = 0;

    // (ix,iy),(kx,ky) are the co-ordinate values of
    // the points
    for (int i = 0; i < poly->nlines; i++) {
        // i and k form a line in polygon
        int ix = poly->lines[i].p1.x, iy = poly->lines[i].p1.y;
        int kx = poly->lines[i].p2.x, ky = poly->lines[i].p2.y;

        // Calculating position of first point
        // w.r.t. clipper line
        int i_pos = (x2 - x1) * (iy - y1) - (y2 - y1) * (ix - x1);

        // Calculating position of second point
        // w.r.t. clipper line
        int k_pos = (x2 - x1) * (ky - y1) - (y2 - y1) * (kx - x1);

        // Case 1 : When both points are inside
        if (i_pos < 0 && k_pos < 0) {
            //Only second point is added
            new_points[new_poly_size][0] = kx;
            new_points[new_poly_size][1] = ky;
            new_poly_size++;
        }

            // Case 2: When only first point is outside
        else if (i_pos >= 0 && k_pos < 0) {
            // Point of intersection with edge
            // and the second point is added
            new_points[new_poly_size][0] = x_intersect(x1,
                                                       y1, x2, y2, ix, iy, kx, ky);
            new_points[new_poly_size][1] = y_intersect(x1,
                                                       y1, x2, y2, ix, iy, kx, ky);
            new_poly_size++;

            new_points[new_poly_size][0] = kx;
            new_points[new_poly_size][1] = ky;
            new_poly_size++;
        }

            // Case 3: When only second point is outside
        else if (i_pos < 0 && k_pos >= 0) {
            //Only point of intersection with edge is added
            new_points[new_poly_size][0] = x_intersect(x1,
                                                       y1, x2, y2, ix, iy, kx, ky);
            new_points[new_poly_size][1] = y_intersect(x1,
                                                       y1, x2, y2, ix, iy, kx, ky);
            new_poly_size++;
        }

            // Case 4: When both points are outside
        else {
            //No points are added
        }
    }
    // Copying new points into original array
    // and changing the no. of vertices
    poly->nlines = new_poly_size;
    for (int i = 0; i < poly->nlines; i++)
    {
        poly->lines[i].p1.x = new_points[i][0];
        poly->lines[i].p1.y = new_points[i][1];

        poly->lines[i].p2.x = new_points[(i+1) % poly->nlines][0];
        poly->lines[i].p2.y = new_points[(i+1) % poly->nlines][1];
    }
}

// Implements Sutherland–Hodgman algorithm
void ClipPolygon(POLYGON* poly, int clipper[][2]) {
    //i and k are two consecutive indexes

    for (int i = 0; i < 4; i++) {
        int k = (i + 1) % 4;

        // We pass the current array of vertices, it's size
        // and the end points of the selected clipper line
        clip(poly,
             clipper[i][0],
             clipper[i][1],
             clipper[k][0],
             clipper[k][1]);
    }
}

void ClipPolygons(int clipper[][2]) {
    for (int i = 0; i < LEN(polygons); ++i) {
        ClipPolygon(polygons[i], clipper);
    }
}

void plot(int x, int y) {
    if (x < 0 || y < 0 || x > H_SIZE - 1 || y > V_SIZE - 1)return;

    int tpos = (y % 255) * 255 + x % 255;

    switch (textureId) {
        case 0:
            buffer[x][y] = global_color;
            break;
        case 2:
            buffer[x][y].r = (float) GT[tpos][0];
            buffer[x][y].g = (float) GT[tpos][1];
            buffer[x][y].b = (float) GT[tpos][2];
            break;
        case 7:
            buffer[x][y].r = (float) PT[tpos][0];
            buffer[x][y].g = (float) PT[tpos][1];
            buffer[x][y].b = (float) PT[tpos][2];
            break;
        case 8:
            buffer[x][y].r = (float) PT[tpos][0];
            buffer[x][y].g = (float) PT[tpos][1];
            buffer[x][y].b = (float) PT[tpos][2];
            break;
        case 4:
            buffer[x][y].r = (float) LT[tpos][0];
            buffer[x][y].g = (float) LT[tpos][1];
            buffer[x][y].b = (float) LT[tpos][2];
            break;
        case 6:
            buffer[x][y].r = (float) AT[tpos][0];
            buffer[x][y].g = (float) AT[tpos][1];
            buffer[x][y].b = (float) AT[tpos][2];
            break;
        case 1:
            buffer[x][y].r = (float) CT[tpos][0];
            buffer[x][y].g = (float) CT[tpos][1];
            buffer[x][y].b = (float) CT[tpos][2];
            break;
        case 3:
            buffer[x][y].r = (float) HT[tpos][0];
            buffer[x][y].g = (float) HT[tpos][1];
            buffer[x][y].b = (float) HT[tpos][2];
            break;
        case 5:
            buffer[x][y].r = (float) SJT[tpos][0];
            buffer[x][y].g = (float) SJT[tpos][1];
            buffer[x][y].b = (float) SJT[tpos][2];
            break;
        default:
            buffer[x][y] = global_color;
            break;

    }

}

void set_color(double r, double g, double b) {
    global_color.r = r;
    global_color.g = g;
    global_color.b = b;
}

LINE setLineValues(LINE line) {
    if (line.p1.y > line.p2.y) {
        int x = line.p1.x;
        int y = line.p1.y;
        line.p1.x = line.p2.x;
        line.p1.y = line.p2.y;
        line.p2.x = x;
        line.p2.y = y;
    }

    //printLine(line);
    if (line.p2.y == line.p1.y) {
        //printLine("y igual");
        line.m = 1;
    } else if (line.p2.x == line.p1.x) {
        //printLine("x igual");
        line.m = 0;
    } else {
        line.m = (float) ((float) ((float) line.p2.y - (float) line.p1.y) /
                          (float) ((float) line.p2.x - (float) line.p1.x));
        //printf("m %f\n", line.m);
    }


    line.b = line.p1.y - line.m * line.p1.x;

    if (line.m != 0) {
        line.delta = -1.0 / (float) line.m;
    } else {
        line.delta = 0;
    }
    //printf("delta %f\n", line.delta);
    return line;
}

POLYGON *NewPolygon() {
    POLYGON *np = (POLYGON *) malloc(sizeof(POLYGON));
    np->nlines = -1;
    np->lines = (LINE *) malloc(sizeof(LINE) * 1000); //puntos
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
        //poly->lines[0] = setLineValues(poly->lines[0]);
        poly->nlines++;
    } else if (poly->nlines == 1) {
        poly->lines[poly->nlines].p1.x = poly->lines[poly->nlines - 1].p2.x;
        poly->lines[poly->nlines].p1.y = poly->lines[poly->nlines - 1].p2.y;
        poly->lines[poly->nlines].p2.x = x;
        poly->lines[poly->nlines].p2.y = y;
        //poly->lines[poly->nlines] = setLineValues(poly->lines[poly->nlines]);
        poly->nlines++;

        poly->lines[poly->nlines].p1.x = poly->lines[poly->nlines - 1].p2.x;
        poly->lines[poly->nlines].p1.y = poly->lines[poly->nlines - 1].p2.y;
        poly->lines[poly->nlines].p2.x = poly->lines[0].p1.x;
        poly->lines[poly->nlines].p2.y = poly->lines[0].p1.y;
        //poly->lines[poly->nlines] = setLineValues(poly->lines[poly->nlines]);
        poly->nlines++;
    } else {
        poly->lines[poly->nlines - 1].p1.x = poly->lines[poly->nlines - 2].p2.x;
        poly->lines[poly->nlines - 1].p1.y = poly->lines[poly->nlines - 2].p2.y;
        poly->lines[poly->nlines - 1].p2.x = x;
        poly->lines[poly->nlines - 1].p2.y = y;
        //poly->lines[poly->nlines-1] = setLineValues(poly->lines[poly->nlines-1]);

        poly->lines[poly->nlines].p1.x = poly->lines[poly->nlines - 1].p2.x;
        poly->lines[poly->nlines].p1.y = poly->lines[poly->nlines - 1].p2.y;
        poly->lines[poly->nlines].p2.x = poly->lines[0].p1.x;
        poly->lines[poly->nlines].p2.y = poly->lines[0].p1.y;
        //poly->lines[poly->nlines] = setLineValues(poly->lines[poly->nlines]);
        poly->nlines++;
    }
}

void printLine(LINE line) {
    printf("p1 %d %d p2 %d %d \n", line.p1.x, line.p1.y, line.p2.x, line.p2.y);
}

void PaintPolygon(POLYGON *poly) {
    if (poly->nlines < 2) { return; }

    LINE auxLines[poly->nlines];

    int initialy = 100000000;
    int finaly = -100000000;
    POINT intersections[50];
    int nIntersections = 0;
    POINT p1, p2;

    int lineState[poly->nlines];
    for (int i = 0; i < poly->nlines; i++) {
        auxLines[i] = setLineValues(poly->lines[i]);

        lineState[i] = 0;
        if (auxLines[i].p1.y < initialy) initialy = auxLines[i].p1.y;
        if (auxLines[i].p2.y > finaly) finaly = auxLines[i].p2.y;
    }

    int vertex = 0;

    for (int i = initialy; i <= finaly; i++) {
        nIntersections = 0;
        vertex = 0;
        for (int j = 0; j < poly->nlines; j++) {
            if (auxLines[j].p1.y != auxLines[j].p2.y) {
                if (auxLines[j].p1.y == i) { //punto arriba
                    //printf("y = %d  una linea activada\n", i);
                    vertex++;
                    lineState[j] = 2;
                } else if (auxLines[j].p2.y == i) { // punto abajo
                    vertex++;
                    lineState[j] = 3;
                } else if (auxLines[j].p1.y < i) {
                    lineState[j] = 1;
                }
                if (auxLines[j].p2.y < i) {
                    //printf("y = %d  una linea desactivada\n", i);
                    lineState[j] = 0;
                }

            } else {
                if (auxLines[j].p1.y == i) {
                    lineState[j] = -1;
                    vertex += 2;

                }
            }
        }


        //arreglar codos y rodillas y gradas

        if (vertex > 1) {
            for (int l = 1; l < poly->nlines; ++l) {
                if (lineState[l - 1] == 2 && lineState[l] == 3) {
                    if (auxLines[l - 1].p1.x == auxLines[l].p2.x) {
                        intersections[nIntersections].x = auxLines[l - 1].p1.x;
                        intersections[nIntersections].y = i;

                        //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                        nIntersections++;

                        lineState[l] = 0;
                        lineState[l - 1] = 0;

                    }
                } else if (lineState[l - 1] == 3 && lineState[l] == 2) {
                    if (auxLines[l - 1].p2.x == auxLines[l].p1.x) {
                        intersections[nIntersections].x = auxLines[l - 1].p2.x;
                        intersections[nIntersections].y = i;

                        //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                        nIntersections++;

                        lineState[l] = 0;
                        lineState[l - 1] = 0;

                    }
                } else if (lineState[l - 1] == 3 && lineState[l] == 3) {
                    if (auxLines[l - 1].p2.x == auxLines[l].p2.x) {
                        intersections[nIntersections].x = auxLines[l - 1].p2.x;
                        intersections[nIntersections].y = i - 1;

                        //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                        nIntersections++;
                        intersections[nIntersections].x = auxLines[l - 1].p1.x;
                        intersections[nIntersections].y = i - 1;

                        //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                        nIntersections++;
                        lineState[l] = 0;
                        lineState[l - 1] = 0;

                    }
                } else if (lineState[l - 1] == 2 && lineState[l] == 2) {
                    if (auxLines[l - 1].p1.x == auxLines[l].p1.x) {
                        intersections[nIntersections].x = auxLines[l - 1].p1.x;
                        intersections[nIntersections].y = i;

                        //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                        nIntersections++;
                        intersections[nIntersections].x = auxLines[l - 1].p1.x;
                        intersections[nIntersections].y = i;

                        //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                        nIntersections++;
                        lineState[l] = 0;
                        lineState[l - 1] = 0;

                    }
                } else if ((lineState[l] == -1 && lineState[l - 1] == 2 && lineState[l + 1] == 2) ||
                           (lineState[l - 1] == -1 && lineState[l - 1] == 3 && lineState[l + 1] == 3)) {
                    lineState[l] = 0;
                    lineState[l - 1] = 0;
                    lineState[l + 1] = 0;
                    //printf("%d\n", j);
                    intersections[nIntersections].x = auxLines[l].p1.x;
                    intersections[nIntersections].y = i;

                    //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                    nIntersections++;
                    //printf("%d\n", j);
                    intersections[nIntersections].x = auxLines[l].p2.x;
                    intersections[nIntersections].y = i;

                    //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                    nIntersections++;
                } else if ((lineState[l] == -1 && lineState[l - 1] == 3 && lineState[l + 1] == 2) ||
                           (lineState[l] == -1 && lineState[l - 1] == 2 && lineState[l + 1] == 3)) {
                    lineState[l] = 0;
                    lineState[l - 1] = 0;
                    lineState[l + 1] = 0;
                    //printf("%d\n", j);
                    intersections[nIntersections].x = auxLines[l].p1.x;
                    intersections[nIntersections].y = i - 1;

                    //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                    nIntersections++;
                    //printf("%d\n", j);
                    intersections[nIntersections].x = auxLines[l].p1.x;
                    intersections[nIntersections].y = i - 1;

                    //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                    nIntersections++;
                    //printf("%d\n", j);
                    intersections[nIntersections].x = auxLines[l].p2.x;
                    intersections[nIntersections].y = i - 1;

                    //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                    nIntersections++;
                    //printf("%d\n", j);
                    intersections[nIntersections].x = auxLines[l].p2.x;
                    intersections[nIntersections].y = i - 1;

                    //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                    nIntersections++;
                }

            }

        }

        for (int j = 0; j < poly->nlines; j++) {

            if (lineState[j] == 1) {
                //printf("%d\n", j);
                intersections[nIntersections].x = round(
                        (float) auxLines[j].p1.x - (float) (i - auxLines[j].p1.y) * (float) auxLines[j].delta);
                intersections[nIntersections].y = i;

                //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                nIntersections++;
            }

        }
        //printf("y = %d  numero de intersecciones %d\n", i , nIntersections );
        qsort(intersections, nIntersections, sizeof(POINT), compare);

        for (int j = 1; j < nIntersections; j += 2) {
            if ((intersections[j - 1].y != i && intersections[j].y != i) &&
                (intersections[j - 1].x == intersections[j].x)) {
                j++;
                if (j < nIntersections) bresenham(intersections[j - 1].x, i, intersections[j].x, i);
            } else if ((intersections[j - 1].y != i && intersections[j].y != i) &&
                       (intersections[j - 1].x != intersections[j].x)) {

                bresenham(intersections[j - 1].x, i, intersections[j].x, i);
                j++;
            } else {
                bresenham(intersections[j - 1].x, i, intersections[j].x, i);
            }


        }

    }


}

void init() {
    textureId = 1;
    lineCount = 0;
    tool = 2;
    drawType = 0;
    Tx = 0;
    Ty = 0;
    alpha = 0;
    Sx = 1;
    Sy = 1;
    polygons[0] = polyCartago = NewPolygon();
    polygons[1] = polyGuanacaste = NewPolygon();
    polygons[2] = polyHeredia = NewPolygon();
    polygons[3] = polyLimon = NewPolygon();
    polygons[4] = polySanJose = NewPolygon();
    polygons[5] = polyAlajuela = NewPolygon();
    polygons[6] = polyPuntarenas = NewPolygon();
    polygons[7] = polyPuntarenas1 = NewPolygon();
    polygons[8] = poly = NewPolygon();
    polygonsAux[0] = NewPolygon();
    polygonsAux[1] = NewPolygon();
    polygonsAux[2] = NewPolygon();
    polygonsAux[3] = NewPolygon();
    polygonsAux[4] = NewPolygon();
    polygonsAux[5] = NewPolygon();
    polygonsAux[6] = NewPolygon();
    polygonsAux[7] = NewPolygon();
    polygonsAux[8] = NewPolygon();
}

void MyKeyboardFunc(unsigned char Key, int x, int y) {
    switch (Key) {
        case '1':
            drawType = 0;
            break;
        case '2':
            drawType = 1;
            break;
        case '3':
            drawType = 2;
            break;
        case '+':
            Sx += 0.05 * r;
            Sy += 0.05 * r;
            break;
        case '-':
            Sx -= 0.05 * r;
            Sy -= 0.05 * r;
            break;
        case 'i':
            alpha += 0.05 * r;
            break;
        case 'j':
            alpha -= 0.05 * r;
            break;
        case 'f':
            r = 3;
            break;
        case 's':
            r = 1;
            break;
    };
    glutPostRedisplay();
}

void SpecialInput(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:

            // Yc += 10;
            Ty += 3 * r;
            break;
        case GLUT_KEY_DOWN:

            //Yc -= 10;
            Ty -= 3 * r;
            break;
        case GLUT_KEY_LEFT:
            //Xc += 10;
            Tx += 3 * r;
            break;
        case GLUT_KEY_RIGHT:
            //Xc -= 10;
            Tx -= 3 * r;
            break;
    };
    glutPostRedisplay();
}

void myMouseFunc(int button, int state, int x, int y) {

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

    if (y0 == y1) {
        Xp = x0;
        Yp = y0;
        plot(Xp, Yp);
        while (Xp != x1) {
            Xp++;
            plot(Xp, Yp);
        }
    } else if (x0 == x1) {

        if (y0 > y1) {
            int aux = y0;
            y0 = y1;
            y1 = aux;
        }
        Xp = x0;
        Yp = y0;
        plot(Xp, Yp);
        while (Yp != y1) {
            Yp++;
            plot(Xp, Yp);
        }
    } else {
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

}

void clear_scene() {
    int i, j, tpos;
    for (i = 0; i < H_SIZE; i++) {
        for (j = 0; j < V_SIZE; j++) {
            tpos = (j % 255) * 255 + i % 255;
            buffer[i][j].r = SEA[tpos][0];
            buffer[i][j].g = SEA[tpos][1];
            buffer[i][j].b = SEA[tpos][2];
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

void DrawPolygon(POLYGON *poly) {
    for (int i = 0; i < poly->nlines; i++) {
        set_color(0, 0, 0);
        bresenham(
                poly->lines[i].p1.x,
                poly->lines[i].p1.y,
                poly->lines[i].p2.x,
                poly->lines[i].p2.y
        );
    }
}

void DrawPolygons() {
    textureId = 0;
    set_color(0, 0, 0);
    for (int i = 0; i < LEN(polygons); ++i) {
        DrawPolygon(polygons[i]);
    }
}

void TexturePolygons() {
    for (int i = 0; i < LEN(polygons);
    ++i) {
        textureId = i + 1;
        PaintPolygon(polygons[i]);
    }
}

void PaintPolygons() {
    textureId = 0;

    set_color(0.5, 0, 0.5);
    PaintPolygon(polygons[0]);

    set_color(0, 0.5, 0.2);
    PaintPolygon(polygons[1]);

    set_color(0.8, 0.3, 0);
    PaintPolygon(polygons[2]);

    set_color(0.1, 0.2, 0.3);
    PaintPolygon(polygons[3]);

    set_color(0.4, 0.3, 0.9);
    PaintPolygon(polygons[4]);

    set_color(0.1, 0.9, 0.1);
    PaintPolygon(polygons[5]);

    set_color(0.9, 0.2, 0.1);
    PaintPolygon(polygons[6]);

    set_color(0.9, 0.2, 0.1);
    PaintPolygon(polygons[7]);

}

void reset_db() {
    for (int i = 0; i < LEN(polygons); ++i) {
        polygons[i]->nlines = 1000;
        for (int j = 0; j < polygons[i]->nlines; ++j) {
            polygons[i]->lines[j].p1.x = polygonsAux[i]->lines[j].p1.x;
            polygons[i]->lines[j].p1.y = polygonsAux[i]->lines[j].p1.y;
            polygons[i]->lines[j].p2.x = polygonsAux[i]->lines[j].p2.x;
            polygons[i]->lines[j].p2.y = polygonsAux[i]->lines[j].p2.y;
        }
    }

}

void renderScene(void) {
    clear_scene();
    reset_db();

    int clipper_points[][2] = {{150,150},
                               {150, H_SIZE - 150},
                               {V_SIZE -150, H_SIZE - 150},
                               {V_SIZE -150, 150} };

    bresenham(clipper_points[0][0], clipper_points[0][1], clipper_points[1][0], clipper_points[1][1]);
    bresenham(clipper_points[1][0], clipper_points[1][1], clipper_points[2][0], clipper_points[2][1]);
    bresenham(clipper_points[2][0], clipper_points[2][1], clipper_points[3][0], clipper_points[3][1]);
    bresenham(clipper_points[3][0], clipper_points[3][1], clipper_points[0][0], clipper_points[0][1]);
    ClipPolygons(clipper_points);

    S_db(Sx, Sy);
    R_db(alpha);
    T_db(Tx, Ty);

    if (drawType == 0)DrawPolygons();
    if (drawType == 1)PaintPolygons();
    if (drawType == 2)TexturePolygons();

    draw_scene();
}

void initTextures() {
    FILE *streamIn;

    int byte, i;

    //Guanacaste

    streamIn = fopen("./GT.tga", "r");
    if (streamIn == (FILE *) 0) {
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    GT = malloc(255 * 255 * sizeof(int *));


    for (i = 0; i < 18; i++) byte = getc(streamIn);
    for (i = 0; i < 255 * 255; i++) {    // foreach pixel
        GT[i] = malloc(3 * sizeof(int));
        GT[i][2] = (float) getc(streamIn) / 255.0;  // use BMP 24bit with no alpha channel
        GT[i][1] = (float) getc(streamIn) / 255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        GT[i][0] = (float) getc(streamIn) / 255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Puntarenas

    streamIn = fopen("./PT.tga", "r");
    if (streamIn == (FILE *) 0) {
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    PT = malloc(255 * 255 * sizeof(int *));


    for (i = 0; i < 18; i++) byte = getc(streamIn);
    for (i = 0; i < 255 * 255; i++) {    // foreach pixel
        PT[i] = malloc(3 * sizeof(int));
        PT[i][2] = (float) getc(streamIn) / 255.0;  // use BMP 24bit with no alpha channel
        PT[i][1] = (float) getc(streamIn) / 255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        PT[i][0] = (float) getc(streamIn) / 255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Limon

    streamIn = fopen("./LT.tga", "r");
    if (streamIn == (FILE *) 0) {
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    LT = malloc(255 * 255 * sizeof(int *));


    for (i = 0; i < 18; i++) byte = getc(streamIn);
    for (i = 0; i < 255 * 255; i++) {    // foreach pixel
        LT[i] = malloc(3 * sizeof(int));
        LT[i][2] = (float) getc(streamIn) / 255.0;  // use BMP 24bit with no alpha channel
        LT[i][1] = (float) getc(streamIn) / 255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        LT[i][0] = (float) getc(streamIn) / 255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Alajuela

    streamIn = fopen("./AT.tga", "r");
    if (streamIn == (FILE *) 0) {
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    AT = malloc(255 * 255 * sizeof(int *));


    for (i = 0; i < 18; i++) byte = getc(streamIn);
    for (i = 0; i < 255 * 255; i++) {    // foreach pixel
        AT[i] = malloc(3 * sizeof(int));
        AT[i][2] = (float) getc(streamIn) / 255.0;  // use BMP 24bit with no alpha channel
        AT[i][1] = (float) getc(streamIn) / 255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        AT[i][0] = (float) getc(streamIn) / 255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Guanacaste

    streamIn = fopen("./SJT.tga", "r");
    if (streamIn == (FILE *) 0) {
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    SJT = malloc(255 * 255 * sizeof(int *));


    for (i = 0; i < 18; i++) byte = getc(streamIn);
    for (i = 0; i < 255 * 255; i++) {    // foreach pixel
        SJT[i] = malloc(3 * sizeof(int));
        SJT[i][2] = (float) getc(streamIn) / 255.0;  // use BMP 24bit with no alpha channel
        SJT[i][1] = (float) getc(streamIn) / 255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        SJT[i][0] = (float) getc(streamIn) / 255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Cartago

    streamIn = fopen("./CT.tga", "r");
    if (streamIn == (FILE *) 0) {
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    CT = malloc(255 * 255 * sizeof(int *));


    for (i = 0; i < 18; i++) byte = getc(streamIn);
    for (i = 0; i < 255 * 255; i++) {    // foreach pixel
        CT[i] = malloc(3 * sizeof(int));
        CT[i][2] = (float) getc(streamIn) / 255.0;  // use BMP 24bit with no alpha channel
        CT[i][1] = (float) getc(streamIn) / 255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        CT[i][0] = (float) getc(streamIn) / 255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Heredia

    streamIn = fopen("./HT.tga", "r");
    if (streamIn == (FILE *) 0) {
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    HT = malloc(255 * 255 * sizeof(int *));


    for (i = 0; i < 18; i++) byte = getc(streamIn);
    for (i = 0; i < 255 * 255; i++) {    // foreach pixel
        HT[i] = malloc(3 * sizeof(int));
        HT[i][2] = (float) getc(streamIn) / 255.0;  // use BMP 24bit with no alpha channel
        HT[i][1] = (float) getc(streamIn) / 255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        HT[i][0] = (float) getc(streamIn) / 255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //SEA

    streamIn = fopen("./SEA.tga", "r");
    if (streamIn == (FILE *) 0) {
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    SEA = malloc(255 * 255 * sizeof(int *));


    for (i = 0; i < 18; i++) byte = getc(streamIn);
    for (i = 0; i < 255 * 255; i++) {    // foreach pixel
        SEA[i] = malloc(3 * sizeof(int));
        SEA[i][2] = (float) getc(streamIn) / 255.0;  // use BMP 24bit with no alpha channel
        SEA[i][1] = (float) getc(streamIn) / 255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        SEA[i][0] = (float) getc(streamIn) / 255.0;  // reverse-order array indexing fixes RGB issue...

    }
}

int main(int argc, char **argv) {
    int i, j, length;


    initTextures();
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

    char *rutas[] = {
            "mapas/cartago.txt",
            "mapas/guanacaste.txt",
            "mapas/heredia.txt",
            "mapas/limon.txt",
            "mapas/sanjose.txt",
            "mapas/alajuela.txt",
            "mapas/puntarenas.txt",
            "mapas/puntarenas1.txt",
    };
    for (int k = 0; k < LEN(rutas);
    ++k) {
        file = fopen(rutas[k], "r");

        int x;
        int y;

        i = 0;
        char line[25];
        while (fgets(line, 25, file)) {
            //if(line == NULL)break;
            // double row[ssParams->nreal + 1];
            char *tmp = strdup(line);

            int j = 0;
            const char *tok;
            for (tok = strtok(line, ","); tok && *tok; j++, tok = strtok(NULL, "\t\n")) {
                if (j == 0) x = atoi(tok);
                else y = atoi(tok);
            }

            AddPolygonLine(polygons[k], x, y);
            AddPolygonLine(polygonsAux[k], x, y);

            free(tmp);
            i++;
        }
        //T_polygon(polygons[k], 100, 100);
        //S_polygon(polygons[k], 1.2, 1.2);

    }

    // enter GLUT event processing cycle
    glutMainLoop();

    printf("Fin del programa %s...\n\n", argv[0]);


    return 1;
}