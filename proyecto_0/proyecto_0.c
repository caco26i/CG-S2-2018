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
#include <stdbool.h>
#include <malloc.h>
#include <GL/glut.h>

#include "proyecto_0.h"

void draw_scene();

int main(int argc, char **argv) {
    int i, j, length;


    int edad;

    if (argc <= 3) {
        printf("Debes ingresar mas parametros...\n");
        return 1;
    }

    printf("Resolución: %s\n", argv[1]);
    printf("# cantidad_lineas: %s\n", argv[2]);
    printf("# cantidad_veces: %s\n", argv[3]);

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
            printf("Ingrese un número de cantidad_lineas válido\n");
            return 1;
        }
    sscanf(argv[2], "%d", &cantidad_lineas);

    length = strlen(argv[3]);
    for (i = 0; i < length; i++)
        if (!isdigit(argv[3][i])) {
            printf("Ingrese un número de cantidad_veces válido\n");
            return 1;
        }
    sscanf(argv[3], "%d", &cantidad_veces);

    buffer = (COLOR **) malloc(H_SIZE * sizeof(COLOR *));
    buffer_random_lines = (LINE *) malloc(cantidad_lineas * sizeof(LINE *));


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
    draw_lines();

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

void set_color(double r, double g, double b) {
    global_color.r = r;
    global_color.g = g;
    global_color.b = b;
}

void plot(int x, int y) {
    buffer[x][y] = global_color;
}

void generate_random_lines() {
    int i;
    for (i = 0; i < cantidad_lineas; i++) {
        buffer_random_lines[i].x0 = 1 + get_int_random(H_SIZE - 2);
        buffer_random_lines[i].x1 = 1 + get_int_random(H_SIZE - 2);
        buffer_random_lines[i].y0 = 1 + get_int_random(H_SIZE - 2);
        buffer_random_lines[i].y1 = 1 + get_int_random(H_SIZE - 2);
    }
}

void line(int x0, int y0, int x1, int y1) {
    long double mx, my, bx, by, x, y;
    int i;
    bool more_x = abs(x1 - x0) > abs(y1 - y0);

    if (more_x) {
        if(x0 > x1) {
            int aux = x0;
            x0 = x1;
            x1 = aux;
            aux = y0;
            y0 = y1;
            y1 = aux;
        }
    } else {
        if(y0 > y1) {
            int aux = x0;
            x0 = x1;
            x1 = aux;
            aux = y0;
            y0 = y1;
            y1 = aux;
        }
    }

    mx = (long double) (y1 - y0) / (x1 - x0);
    my = (long double) (x1 - x0) / (y1 - y0);
    bx = y0 - mx * x0;
    by = x0 - my * y0;

    if (more_x)
        for (i = x0; i <= x1; i++) {
            y = mx * i + bx;
            plot(i, round(y));
        }
    else
        for (i = y0; i <= y1; i++) {
            x = my * i + by;
            plot(round(x), i);
        }

}

void line2(int x0, int y0, int x1, int y1) {
    long double mx, my, y, x;
    int i;

    bool more_x = abs(x1 - x0) > abs(y1 - y0);

    if (more_x) {
        if(x0 > x1) {
            int aux = x0;
            x0 = x1;
            x1 = aux;
            aux = y0;
            y0 = y1;
            y1 = aux;
        }
    } else {
        if(y0 > y1) {
            int aux = x0;
            x0 = x1;
            x1 = aux;
            aux = y0;
            y0 = y1;
            y1 = aux;
        }
    }

    mx = (long double) (y1 - y0) / (x1 - x0);
    my = (long double) (x1 - x0) / (y1 - y0);
    y = y0;
    x = x0;

    if (more_x)
        for (i = x0; i <= x1; i++) {
            plot(i, round(y));
            y += mx;
        }
    else
        for (i = y0; i <= y1; i++) {
            plot(round(x), i);
            x += my;
        }
}

void line3(int x0, int y0, int x1, int y1) {
    long double x, y, paso_x, paso_y;
    int i, ancho;

    ancho = max(abs(x1 - x0), abs(y1 - y0));

    paso_x = (long double) (x1 - x0) / ancho;
    paso_y = (long double) (y1 - y0) / ancho;

    x = x0;
    y = y0;

    for (i = 0; i <= ancho; i++) {
        plot(round(x), round(y));
        x += paso_x;
        y += paso_y;
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
    Delta_E = 2 * dify;
    Delta_NE = 2 * dify - 2 * difx;
    Delta_N = -2 * difx;
    Delta_SE = 2 * difx + 2 * dify;
    Delta_S = 2 * difx;
    Xp = x0;
    Yp = y0;

    if (y0 <= y1) {
        if (abs(difx) >= abs(dify)) {   //Cuadrante 1
            Delta_AX = 1;
            Delta_AY = 0;
            Delta_BX = 1;
            Delta_BY = 1;
            Delta_A = Delta_E;
            Delta_B = Delta_NE;
            d = 2 * dify - difx;
        } else {                        //Cuadrante 2
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
            Delta_BX = 1;
            Delta_BY = 0;
            Delta_AX = 1;
            Delta_AY = -1;
            Delta_B = Delta_E;
            Delta_A = Delta_SE;
            d = 2 * dify + difx;
        } else {                        //Cuadrante 7
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

    while (Xp < x1 || Yp < y1) {
        if (d <= 0) {
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

void draw_lines() {
    int i, vez;
    clock_t tiempo_inicio, tiempo_final;
    double segundos;

    set_color(0, 0, 1);
    tiempo_inicio = clock();
    for (vez = 0; vez < cantidad_veces; vez++)
        for (i = 0; i < cantidad_lineas; i++)
            line(
                    buffer_random_lines[i].x0,
                    buffer_random_lines[i].x1,
                    buffer_random_lines[i].y0,
                    buffer_random_lines[i].y1
            );
    tiempo_final = clock();
    segundos = (double) (tiempo_final - tiempo_inicio) / CLOCKS_PER_SEC;
    printf("algoritmo 1 %f segundos\n\n", segundos);

    set_color(1, 0, 0);
    tiempo_inicio = clock();
    for (vez = 0; vez < cantidad_veces; vez++)
        for (i = 0; i < cantidad_lineas; i++)
            line2(
                    buffer_random_lines[i].x0,
                    buffer_random_lines[i].x1,
                    buffer_random_lines[i].y0,
                    buffer_random_lines[i].y1
            );
    tiempo_final = clock();
    segundos = (double) (tiempo_final - tiempo_inicio) / CLOCKS_PER_SEC;
    printf("algoritmo 2 %f segundos\n\n", segundos);


    set_color(0, 1, 0);
    tiempo_inicio = clock();
    for (vez = 0; vez < cantidad_veces; vez++)
        for (i = 0; i < cantidad_lineas; i++)
            line3(
                    buffer_random_lines[i].x0,
                    buffer_random_lines[i].x1,
                    buffer_random_lines[i].y0,
                    buffer_random_lines[i].y1
            );
    tiempo_final = clock();
    segundos = (double) (tiempo_final - tiempo_inicio) / CLOCKS_PER_SEC;
    printf("algoritmo 3 %f segundos\n\n", segundos);

    set_color(1, 1, 1);
    tiempo_inicio = clock();
    for (vez = 0; vez < cantidad_veces; vez++)
        for (i = 0; i < cantidad_lineas; i++)
            bresenham(
                    buffer_random_lines[i].x0,
                    buffer_random_lines[i].x1,
                    buffer_random_lines[i].y0,
                    buffer_random_lines[i].y1
            );
    tiempo_final = clock();
    segundos = (double) (tiempo_final - tiempo_inicio) / CLOCKS_PER_SEC;
    printf("algoritmo bresenham %f segundos\n\n", segundos);
}



