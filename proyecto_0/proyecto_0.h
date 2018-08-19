/*
 * Instituto Tecnologico de Costa Rica
 * Escuela de Ingenieria en Computacion
 * Computer Graphics
 *
 * Programa: Mesa Example
 * Archivo:  mesa_example.h
 */

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

#define H_SIZE_default 400
#define V_SIZE_default 400

typedef struct {
    double r;
    double g;
    double b;
} COLOR;

int H_SIZE, V_SIZE;
int lineas, veces;

COLOR **buffer;

double r, g, b;

int get_int_random(int);

void clear_scene(void);

void set_color(double, double, double);

void plot(int, int);

void generate_random_lines(void);

void line(int, int, int, int);

void line2(int, int, int, int);

void line3(int, int, int, int);




