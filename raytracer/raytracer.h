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
    POINT point;
    void* object;
    bool set;
} INTERSECTION;

POINT eye;
COLOR **buffer;
COLOR color;
COLOR background;

int window;
int H_SIZE, V_SIZE;


void draw_scene();
void renderScene();
