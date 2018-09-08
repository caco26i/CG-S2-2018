/*
 * Instituto Tecnologico de Costa Rica
 * Escuela de Ingenieria en Computacion
 * Computer Graphics
 *
 * Programa: Proyecto 1
 * Archivo:  proyecto_1.h
 */

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#define H_SIZE_default 400
#define V_SIZE_default 400

typedef struct {
    double r;
    double g;
    double b;
} COLOR;

typedef struct {
    int x;
    int y;
} POINT;

typedef struct {
    POINT p1;
    POINT p2;
    int m;
    int b;
    int delta;
} LINE;

typedef struct {
    int nlines;
    LINE *lines;
} POLYGON;

void bresenham(int, int, int, int);

LINE lines[1000];
int lineCount;
int tool;
POLYGON *poly;
