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

#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))


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
    float m;
    float b;
    float delta;
} LINE;

typedef struct {
    int nlines;
    LINE *lines;
} POLYGON;

typedef struct {
    int npolygons;
    POLYGON *polygon;
} POLYGONS;

COLOR **buffer;
COLOR global_color;

int H_SIZE, V_SIZE;

LINE lines[1000];
int lineCount;
int tool;

POLYGON *poly;
POLYGON *polySanJose;
POLYGON *polyHeredia;
POLYGON *polyAlajuela;
POLYGON *polyGuanacaste;
POLYGON *polyCartago;
POLYGON *polyPuntarenas;
POLYGON *polyPuntarenas1;
POLYGON *polyLimon;
POLYGON *polygons[9];
POLYGON *polygonsAux[9];

long double R_matrix[3][3];

void bresenham(int, int, int, int);
void draw_scene();
void renderScene();
LINE setLineValues(LINE);
