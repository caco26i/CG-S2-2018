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

typedef struct {
    int x0;
    int x1;
    int y0;
    int y1;
} LINE;

int H_SIZE, V_SIZE;
int cantidad_lineas, cantidad_veces;

COLOR **buffer;

LINE *buffer_random_lines;

COLOR global_color;

int get_int_random(int);

void clear_scene(void);

void set_color(double, double, double);

void plot(int, int);

void generate_random_lines(void);

void line(int, int, int, int);

void line2(int, int, int, int);

void line3(int, int, int, int);

void bresenham(int, int, int, int);

void draw_lines();




