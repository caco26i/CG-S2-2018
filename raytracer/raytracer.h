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
#define T_SPHERE 1
#define N_SPHERES 1
#define INF 10000000

typedef struct {
    float x;
    float y;
    float z;
} POINT;

typedef struct {
    double R;
    double G;
    double B;
} COLOR;

typedef struct {
    int type;
    float radius;
    POINT center;
    COLOR color;
} SPHERE;

typedef struct {
    int t;
    void* object;
    bool set;
} INTERSECTION;

typedef struct {
    POINT pmin;
    POINT pmax;
} VIEWPORT;

POINT eye;
COLOR **buffer;
COLOR color;
COLOR background;
SPHERE **spheres;
VIEWPORT viewport;

int window;
int H_SIZE, V_SIZE;


void draw_scene();
void renderScene();
