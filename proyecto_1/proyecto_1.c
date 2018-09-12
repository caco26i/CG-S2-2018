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
#include <IL/il.h>
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

float** GT;
float** PT;
float** LT;
float** AT;
float** CT;
float** SJT;
float** HT;
float** SEA;
int textureId;
int drawType;
int compare (const void * a, const void * b)
{

  POINT *orderA = (POINT *)a;
  POINT *orderB = (POINT *)b;

  return ( orderB->x - orderA->x );
}


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

void R_db(int alpha){
    for (int i = 0; i <  LEN(polygons); ++i) {
        R_polygon(polygons[i],alpha);
    }
}
void T_db(long double Dx, long double Dy){
    for (int i = 0; i <  LEN(polygons); ++i) {
        T_polygon(polygons[i],Dx, Dy);
    }
}
void S_db(long double Sx, long double Sy){
    for (int i = 0; i <  LEN(polygons); ++i) {
        S_polygon(polygons[i],Sx, Sy);
    }
}

void plot(int x, int y) {


    int tpos = (y%255)*255 + x%255;

    switch(textureId){
        case 0:
            buffer[x][y] = global_color;
            break;
        case 2:
            buffer[x][y].r = (float)GT[tpos][0];
            buffer[x][y].g = (float)GT[tpos][1];
            buffer[x][y].b = (float)GT[tpos][2];
            break;
        case 7:
            buffer[x][y].r = (float)PT[tpos][0];
            buffer[x][y].g = (float)PT[tpos][1];
            buffer[x][y].b = (float)PT[tpos][2];
            break;
        case 8:
            buffer[x][y].r = (float)PT[tpos][0];
            buffer[x][y].g = (float)PT[tpos][1];
            buffer[x][y].b = (float)PT[tpos][2];
            break;
        case 4:
            buffer[x][y].r = (float)LT[tpos][0];
            buffer[x][y].g = (float)LT[tpos][1];
            buffer[x][y].b = (float)LT[tpos][2];
            break;
        case 6:
            buffer[x][y].r = (float)AT[tpos][0];
            buffer[x][y].g = (float)AT[tpos][1];
            buffer[x][y].b = (float)AT[tpos][2];
            break;
        case 1:
            buffer[x][y].r = (float)CT[tpos][0];
            buffer[x][y].g = (float)CT[tpos][1];
            buffer[x][y].b = (float)CT[tpos][2];
            break;
        case 3:
            buffer[x][y].r = (float)HT[tpos][0];
            buffer[x][y].g = (float)HT[tpos][1];
            buffer[x][y].b = (float)HT[tpos][2];
            break;
        case 5:
            buffer[x][y].r = (float)SJT[tpos][0];
            buffer[x][y].g = (float)SJT[tpos][1];
            buffer[x][y].b = (float)SJT[tpos][2];
            break;

    }

}

void set_color(double r, double g, double b) {
    global_color.r = r;
    global_color.g = g;
    global_color.b = b;
}

LINE setLineValues(LINE line) {
    if(line.p1.y > line.p2.y){
        int x = line.p1.x;
        int y = line.p1.y;
        line.p1.x = line.p2.x;
        line.p1.y = line.p2.y;
        line.p2.x = x;
        line.p2.y = y;
    }

    //printLine(line);
    if(line.p2.y == line.p1.y){
        //printLine("y igual");
        line.m = 0;
    }else if (line.p2.x == line.p1.x){
        //printLine("x igual");
        line.m = 1;
    }
    else{
        line.m = (float)((float)((float)line.p2.y - (float)line.p1.y)/(float)((float)line.p2.x - (float)line.p1.x));
        //printf("m %f\n", line.m);
    }


    line.b = line.p1.y - line.m * line.p1.x;

    if (line.m != 0) {
        line.delta = -1.0 / (float)line.m;
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
    }else {
        poly->lines[poly->nlines-1].p1.x = poly->lines[poly->nlines - 2].p2.x;
        poly->lines[poly->nlines-1].p1.y = poly->lines[poly->nlines - 2].p2.y;
        poly->lines[poly->nlines-1].p2.x = x;
        poly->lines[poly->nlines-1].p2.y = y;
        //poly->lines[poly->nlines-1] = setLineValues(poly->lines[poly->nlines-1]);

        poly->lines[poly->nlines].p1.x = poly->lines[poly->nlines - 1].p2.x;
        poly->lines[poly->nlines].p1.y = poly->lines[poly->nlines - 1].p2.y;
        poly->lines[poly->nlines].p2.x = poly->lines[0].p1.x;
        poly->lines[poly->nlines].p2.y = poly->lines[0].p1.y;
        //poly->lines[poly->nlines] = setLineValues(poly->lines[poly->nlines]);
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
}
void printLine(LINE line){
    printf("p1 %d %d p2 %d %d \n", line.p1.x,line.p1.y,line.p2.x,line.p2.y);
}

void PaintPolygon(POLYGON * poly){
	if(poly->nlines < 2){return;}

    LINE auxLines[poly->nlines];

	int initialy = 100000000;
	int finaly = -100000000;
	POINT intersections[50];
    int nIntersections = 0;
	POINT p1,p2;

	int lineState[poly->nlines];
	for(int i = 0; i < poly->nlines; i++){
        auxLines[i] = setLineValues(poly->lines[i]);

		lineState[i] = 0;
		if(auxLines[i].p1.y < initialy) initialy = auxLines[i].p1.y;
		if(auxLines[i].p2.y > finaly) finaly = auxLines[i].p2.y;
	}



	for(int i = initialy ; i <= finaly; i++){
        nIntersections = 0;

        for(int j = 0; j < poly->nlines; j++){
            if(auxLines[j].p1.y != auxLines[j].p2.y){
                if(auxLines[j].p1.y == i){
                    //printf("y = %d  una linea activada\n", i);
                    lineState[j] = 2;
                }else if(auxLines[j].p2.y == i){
                    lineState[j] = 3;
                }else if (auxLines[j].p1.y < i) {
                    lineState[j] = 1;
                }
                if(auxLines[j].p2.y < i){
                    //printf("y = %d  una linea desactivada\n", i);
                    lineState[j] = 0;
                }

            }
        }

        int l = 0;
        int k = 0;

        while(l < poly->nlines){
            if(lineState[l] > 1){
                k = l+1;
                while(k < poly->nlines){
                    if(lineState[k] > 1 && l != k){

                        if(lineState[l] == 2 && lineState[k] == 3){
                            if(auxLines[l].p1.x == auxLines[k].p2.x){
                                lineState[k] = 0;
                            }
                        }else if(lineState[l] == 3 && lineState[k] == 2){
                            if(auxLines[l].p2.x == auxLines[k].p1.x){
                                lineState[k] = 0;
                            }
                        }
                    }
                    k++;
                }
            }
            l++;
        }

		for(int j = 0; j < poly->nlines; j++){

            if(lineState[j] > 0){
                //printf("%d\n", j);
                intersections[nIntersections].x = round(auxLines[j].p1.x - (i - auxLines[j].p1.y) * auxLines[j].delta);
                intersections[nIntersections].y = i;

                //printf("x %d y %d\n",intersections[nIntersections].x,intersections[nIntersections].y);
                nIntersections++;
            }

        }
        //printf("y = %d  numero de intersecciones %d\n", i , nIntersections );
        qsort( intersections, nIntersections, sizeof(POINT), compare );

        for (int j = 1; j < nIntersections; j+=2)
        {

            bresenham(intersections[j-1].x,intersections[j-1].y,intersections[j].x,intersections[j].y);

        }

	}


}


void init() {
    textureId = 1;
    lineCount = 0;
    tool = 2;
    drawType = 0;
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
            drawType = 0;
            break;
        case '2':
            drawType = 1;
            break;
        case '3':
            drawType = 2;
            break;
        case '+':
            S_db(1.1, 1.1);
            break;
        case '-':
            S_db(0.9, 0.9);
            break;
        case 'i':
            R_db(1);
            break;
        case 'j':
            R_db(-1);
            break;
    };
    glutPostRedisplay();
}


void SpecialInput(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:
            T_db(0, -10);
            break;
        case GLUT_KEY_DOWN:
            T_db(0, 10);
            break;
        case GLUT_KEY_LEFT:
            T_db(-10, 0);
            break;
        case GLUT_KEY_RIGHT:
            T_db(10, 0);
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
    int i, j,tpos;
    for (i = 0; i < H_SIZE; i++) {
        for (j = 0; j < V_SIZE; j++) {
            tpos = (j%255)*255 + i%255;
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

void DrawPolygons(){
    textureId = 0;
    set_color(0,0,0);
    for (int i = 0; i < LEN(polygons); ++i) {

        DrawPolygon(polygons[i]);
    }
}

void TexturePolygons(){
    for (int i = 0; i < LEN(polygons); ++i) {
        textureId = i+1;
        PaintPolygon(polygons[i]);
    }
}
void PaintPolygons(){
    textureId = 0;

    set_color(0.5,0,0.5);
    PaintPolygon(polygons[0]);

    set_color(0,0.5,0.2);
    PaintPolygon(polygons[1]);

    set_color(0.8,0.3,0);
    PaintPolygon(polygons[2]);
    
    set_color(0.1,0.2,0.3);
    PaintPolygon(polygons[3]);
    
    set_color(0.4,0.3,0.9);
    PaintPolygon(polygons[4]);
    
    set_color(0.1,0.9,0.1);
    PaintPolygon(polygons[5]);

    set_color(0.9,0.2,0.1);
    PaintPolygon(polygons[6]);

    set_color(0.9,0.2,0.1);
    PaintPolygon(polygons[7]);
    
}
void renderScene(void) {
    clear_scene();
    if(drawType == 0)DrawPolygons();
    if(drawType == 1)PaintPolygons();
    if(drawType == 2)TexturePolygons();

    draw_scene();
}

void initTextures(){
    FILE *streamIn;

    int byte,i;

    //Guanacaste

    streamIn = fopen("./GT.tga", "r");
    if (streamIn == (FILE *)0){
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    GT = malloc(255*255 * sizeof(int*));




    for(i=0;i<18;i++) byte = getc(streamIn);
    for(i=0;i<255*255;i++){    // foreach pixel
        GT[i] = malloc(3*sizeof(int));
        GT[i][2] = (float)getc(streamIn)/255.0;  // use BMP 24bit with no alpha channel
        GT[i][1] = (float)getc(streamIn)/255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        GT[i][0] = (float)getc(streamIn)/255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Puntarenas

    streamIn = fopen("./PT.tga", "r");
    if (streamIn == (FILE *)0){
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    PT = malloc(255*255 * sizeof(int*));




    for(i=0;i<18;i++) byte = getc(streamIn);
    for(i=0;i<255*255;i++){    // foreach pixel
        PT[i] = malloc(3*sizeof(int));
        PT[i][2] = (float)getc(streamIn)/255.0;  // use BMP 24bit with no alpha channel
        PT[i][1] = (float)getc(streamIn)/255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        PT[i][0] = (float)getc(streamIn)/255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Limon

    streamIn = fopen("./LT.tga", "r");
    if (streamIn == (FILE *)0){
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    LT = malloc(255*255 * sizeof(int*));




    for(i=0;i<18;i++) byte = getc(streamIn);
    for(i=0;i<255*255;i++){    // foreach pixel
        LT[i] = malloc(3*sizeof(int));
        LT[i][2] = (float)getc(streamIn)/255.0;  // use BMP 24bit with no alpha channel
        LT[i][1] = (float)getc(streamIn)/255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        LT[i][0] = (float)getc(streamIn)/255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Alajuela

    streamIn = fopen("./AT.tga", "r");
    if (streamIn == (FILE *)0){
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    AT = malloc(255*255 * sizeof(int*));




    for(i=0;i<18;i++) byte = getc(streamIn);
    for(i=0;i<255*255;i++){    // foreach pixel
        AT[i] = malloc(3*sizeof(int));
        AT[i][2] = (float)getc(streamIn)/255.0;  // use BMP 24bit with no alpha channel
        AT[i][1] = (float)getc(streamIn)/255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        AT[i][0] = (float)getc(streamIn)/255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Guanacaste

    streamIn = fopen("./SJT.tga", "r");
    if (streamIn == (FILE *)0){
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    SJT = malloc(255*255 * sizeof(int*));




    for(i=0;i<18;i++) byte = getc(streamIn);
    for(i=0;i<255*255;i++){    // foreach pixel
        SJT[i] = malloc(3*sizeof(int));
        SJT[i][2] = (float)getc(streamIn)/255.0;  // use BMP 24bit with no alpha channel
        SJT[i][1] = (float)getc(streamIn)/255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        SJT[i][0] = (float)getc(streamIn)/255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Cartago

    streamIn = fopen("./CT.tga", "r");
    if (streamIn == (FILE *)0){
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    CT = malloc(255*255 * sizeof(int*));




    for(i=0;i<18;i++) byte = getc(streamIn);
    for(i=0;i<255*255;i++){    // foreach pixel
        CT[i] = malloc(3*sizeof(int));
        CT[i][2] = (float)getc(streamIn)/255.0;  // use BMP 24bit with no alpha channel
        CT[i][1] = (float)getc(streamIn)/255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        CT[i][0] = (float)getc(streamIn)/255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //Heredia

    streamIn = fopen("./HT.tga", "r");
    if (streamIn == (FILE *)0){
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    HT = malloc(255*255 * sizeof(int*));




    for(i=0;i<18;i++) byte = getc(streamIn);
    for(i=0;i<255*255;i++){    // foreach pixel
        HT[i] = malloc(3*sizeof(int));
        HT[i][2] = (float)getc(streamIn)/255.0;  // use BMP 24bit with no alpha channel
        HT[i][1] = (float)getc(streamIn)/255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        HT[i][0] = (float)getc(streamIn)/255.0;  // reverse-order array indexing fixes RGB issue...

    }

    //SEA

    streamIn = fopen("./SEA.tga", "r");
    if (streamIn == (FILE *)0){
        printf("File opening error ocurred. Exiting program.\n");
        exit(0);
    }

    SEA = malloc(255*255 * sizeof(int*));




    for(i=0;i<18;i++) byte = getc(streamIn);
    for(i=0;i<255*255;i++){    // foreach pixel
        SEA[i] = malloc(3*sizeof(int));
        SEA[i][2] = (float)getc(streamIn)/255.0;  // use BMP 24bit with no alpha channel
        SEA[i][1] = (float)getc(streamIn)/255.0;  // BMP uses BGR but we want RGB, grab byte-by-byte
        SEA[i][0] = (float)getc(streamIn)/255.0;  // reverse-order array indexing fixes RGB issue...

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
            //if(line == NULL)break;
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
                //printf("%d - %d,%d", k, x, y);
                //printf("\n");
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
}
