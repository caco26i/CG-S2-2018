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
#include "raytracer.h"

void plot(int x, int y, COLOR c) {
    if (x < 0 || y < 0 || x > H_SIZE - 1 || y > V_SIZE - 1)return;
    buffer[x][y] = c;
}

void set_color(double r, double g, double b) {
    color.R = r;
    color.G = g;
    color.B = b;
}

void init() {

}

void MyKeyboardFunc(unsigned char Key, int x, int y) {
    switch (Key) {
        case 27: // Escape key
            glutDestroyWindow(window);
            exit(1);
            break;
    };
}

void draw_scene() {
    int i, j;

    for (i = 0; i < H_SIZE; i++) {
        for (j = 0; j < V_SIZE; j++) {
            glColor3f(buffer[i][j].R, buffer[i][j].G, buffer[i][j].B);
            glBegin(GL_POINTS);
            glVertex2i(i, j);
            glEnd();
        }
    }
    glFlush();
}

void renderScene(void) {
    draw_scene();
}

INTERSECTION *First_Intersection(POINT a, POINT d) {
    INTERSECTION *interseccion;
    long double tmin;
    interseccion = NULL;
    tmin = 0;
//    ∀ objeto en la escena
//            {
//                    calcular intersección entre rayo y objeto;
//                    Si hay interseccion y distancia al objeto < tmin
//                    {
//                        tmin = distancia
//                        a
//                        intersección;
//                        interseccion = interseccion
//                        con objeto;
//                    }
//            }
    return interseccion;
}


COLOR De_que_color(POINT a, POINT d) {
    COLOR color;
    INTERSECTION* intersection;
    intersection = First_Intersection(a, d);

    if (!intersection)
        color = background;
    else {
        SPHERE *obj = (SPHERE *) intersection->object;
        color = obj->color;
    }

    return color;
}

void raytracer() {
    int x_min;
    int y_min;
    int x_max;
    int y_max;
    POINT w;
    POINT d;
    for (int i = 0; i < H_SIZE; i++) {
        for (int j = 0; j < V_SIZE; j++) {
            w.x = (long double) (i + 1 / 2) * (x_max - x_min) / H_SIZE + x_min;
            w.y = (long double) (j + 1 / 2) * (y_max - y_min) / V_SIZE + y_min;
            w.z = 0;

            long double L = (long double) sqrt(pow(w.x - eye.x, 2) + pow(w.y - eye.y, 2) + pow(w.z - eye.z, 2));
            d.x = (long double) (w.x - eye.x) / L;
            d.y = (long double) (w.y - eye.y) / L;
            d.z = (long double) (w.z - eye.z) / L;
            COLOR color = De_que_color(eye, d);
            plot(i, j, color);
        }
    }
}

int main(int argc, char **argv) {
    int i, j, length;

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

    set_color(1, 1, 1);

    init();

    glutInit(&argc, argv);

    // init GLUT and create Window
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(H_SIZE, V_SIZE);
    window = glutCreateWindow("Raytracer");
    glutKeyboardFunc(MyKeyboardFunc);

    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluOrtho2D(0.0, H_SIZE, V_SIZE, 0.0);
    // register callbacks
    glutDisplayFunc(renderScene);

/*
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
*/
    // enter GLUT event processing cycle
    glutMainLoop();

    printf("Fin del programa %s...\n\n", argv[0]);


    return 1;
}