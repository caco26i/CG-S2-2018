/*
 * Instituto Tecnologico de Costa Rica
 * Escuela de Ingenieria en Computacion
 * Computer Graphics
 *
 * Programa: Proyecto 0
 * Archivo:  proyecto_0.c
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include "proyecto_0.h"
#include "malloc.h"
#include <GL/glut.h>

void draw_scene();

int main(int argc, char **argv) {
    int i, j, length;


    int edad;

    if (argc <= 3) {
        printf("Debes ingresar mas parametros...\n");
        return 1;
    }

    printf("Resolución: %s\n", argv[1]);
    printf("# lineas: %s\n", argv[2]);
    printf("# veces: %s\n", argv[3]);

    length = strlen(argv[1]);
    for (i = 0; i < length; i++)
        if (!isdigit(argv[1][i])) {
            printf("Ingrese una resulución válida\n");
            return 1;
        }

    // Set Res
    sscanf(argv[1], "%d", &H_SIZE);
    sscanf(argv[1], "%d", &V_SIZE);

    length = strlen(argv[2]);
    for (i = 0; i < length; i++)
        if (!isdigit(argv[2][i])) {
            printf("Ingrese un número de lineas válido\n");
            return 1;
        }
    sscanf(argv[2], "%d", &lineas);

    length = strlen(argv[3]);
    for (i = 0; i < length; i++)
        if (!isdigit(argv[3][i])) {
            printf("Ingrese un número de veces válido\n");
            return 1;
        }
    sscanf(argv[3], "%d", &veces);

    buffer = (COLOR **) malloc(H_SIZE * sizeof(COLOR *));

    for (i = 0; i < H_SIZE; i++) {
        buffer[i] = (COLOR *) malloc(V_SIZE * sizeof(COLOR));
    }

    // Random Seed
    // should only be called once
    time_t t;
    srand((unsigned) time(&t));

    clear_scene();
    set_color(1, 1, 1);
    generate_random_lines();

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(H_SIZE, V_SIZE);
    glutCreateWindow("Proyecto 0");
    glClear(GL_COLOR_BUFFER_BIT);
    gluOrtho2D(-0.5, H_SIZE + 0.5, -0.5, V_SIZE + 0.5);
    glutDisplayFunc(draw_scene);
    glutMainLoop();

    printf("Fin del programa %s...\n\n", argv[0]);
    return 0;
}

int get_int_random(int max_number) {
    return rand() % max_number;
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
    static int last_x = 0;
    int i, j;
    COLOR color;

    for (i = 0; i < last_x; i++) {
        for (j = 0; j < V_SIZE; j++) {
            glColor3f(buffer[i][j].r, buffer[i][j].g, buffer[i][j].b);
            glBegin(GL_POINTS);
            glVertex2i(i, j);
            glEnd();
        }
    }

    for (i = last_x; i < H_SIZE; i++) {
        for (j = 0; j < V_SIZE; j++) {
//            buffer[i][j].r = (double) (i % (H_SIZE / 10)) / (double) (H_SIZE / 10);
//            buffer[i][j].g = (double) (j % (V_SIZE / 10)) / (double) (V_SIZE / 10);
//            buffer[i][j].b = (double) (i) / (double) (H_SIZE);
            glColor3f(buffer[i][j].r, buffer[i][j].g, buffer[i][j].b);
            glBegin(GL_POINTS);
            glVertex2i(i, j);
            glEnd();
            last_x = i;
        }
    }

    glFlush();
}

void set_color(double r_p, double g_p, double b_p) {
    r = r_p;
    g = g_p;
    b = b_p;
}

void plot(int x, int y) {
    buffer[x][y].r = r;
    buffer[x][y].g = g;
    buffer[x][y].b = b;
}

void generate_random_lines() {
    int i;
    for (i = 0; i < lineas; i++)
        line3(get_int_random(H_SIZE), get_int_random(H_SIZE), get_int_random(H_SIZE), get_int_random(H_SIZE));
}

void line(int x0, int y0, int x1, int y1) {
    long double m, b, y;
    int i;

    m = (long double) (y1 - y0) / (x1 - x0);
    b = (long double) y0 - m * x0;

    for (i = x0; i <= x1; i++) {
        y = m * i + b;
        plot(i, round(y));
    }
}

void line2(int x0, int y0, int x1, int y1) {
    long double m, y;
    int i;

    m = (long double) (y1 - y0) / (x1 - x0);
    y = y0;

    for (i = x0; i <= x1; i++) {
        plot(i, round(y));
        y += m;
    }
}

void line3(int x0, int y0, int x1, int y1) {
    long double x, y, paso_x, paso_y;
    int i, ancho;

    ancho = max(abs(x1 - x0), abs(y1 - y0));

    paso_x = (long double) (x1-x0) / ancho;
    paso_y = (long double) (y1-y0) / ancho;

    x = x0;
    y = y0;

    for (i = 0; i <= ancho; i++) {
        plot(round(x), round(y));
        x += paso_x;
        y += paso_y;
    }
}



