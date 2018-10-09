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


void initSpheres(){
    for (int i = 0; i < N_SPHERES; ++i)
    {
        spheres[i]->type = T_SPHERE;
        spheres[i]->radius = 100;
        spheres[i]->center.x = 500;
        spheres[i]->center.y = 500;
        spheres[i]->center.z = 50;
        spheres[i]->color = color;
    }
}

void init() {
    eye.x = 500;
    eye.y = 500;
    eye.z = -10;

    color.R = 0;
    color.G = 0;
    color.B = 0;

    background.R = 1;
    background.G = 0;
    background.B = 0;

    spheres = (SPHERE**) malloc(sizeof(SPHERE*)*N_SPHERES);
    for (int i = 0; i < N_SPHERES; ++i)
    {
        spheres[i] = (SPHERE*)malloc(sizeof(SPHERE));
    }
    initSpheres();

    viewport.pmin.x = 0;
    viewport.pmin.y = 0;
    viewport.pmin.z = 0;
    viewport.pmax.x = 1000;
    viewport.pmax.y = 1000;
    viewport.pmax.z = 0;
}


void MyKeyboardFunc(unsigned char Key, int x, int y) {
    switch (Key) {
        case 27: // Escape key
            glutDestroyWindow(window);
            exit(1);
            break;
    };
}



INTERSECTION IntersectionSphere(SPHERE *sphere, POINT e, POINT d){
    INTERSECTION intersection;
    float a = pow(d.x,2) + pow(d.y,2) + pow(d.z,2);
    float b = 2 * ((e.x - sphere->center.x) + (e.y - sphere->center.y) + (e.z - sphere->center.z));
    float g = pow((e.x - sphere->center.x),2) + pow((e.y - sphere->center.y),2) + pow((e.z - sphere->center.z),2) - pow(sphere->radius,2);
    float delta = pow(b,2) - 4 * g * a;
    if(delta > 0.0){

    }else{
        
    }
    
    if(delta < -0.001){
        intersection.t = INF;
    }else if(delta < 0.001){
        intersection.t = -b/2*a;

    }else{
        float t1,t2;
        t1 = (-b + sqrt(delta))/2*a;
        t2 = (-b - sqrt(delta))/2*a;
        if(t1 < -0.001 && t2 < -0.001){
            intersection.t = INF;
        }else if(t1 < -0.001){
            intersection.t = t2;
        }else if(t2 < -0.001){
            intersection.t = t1;
        }else{
            intersection.t = (t1<t2)?t1:t2;
        }
    }
    if(intersection.t < -0.001)intersection.t = INF;
    intersection.object = (void*)sphere;
    return intersection;
}

INTERSECTION First_Intersection(POINT e, POINT d) {
    long double tmin;
    INTERSECTION intersection;
    intersection.t = INF;
    INTERSECTION auxIntersection;
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

    for (int i = 0; i < N_SPHERES; ++i)
    {
        auxIntersection = IntersectionSphere(spheres[i], e, d);
        if(auxIntersection.t >= 0 && auxIntersection.t < intersection.t)intersection = auxIntersection;
    }
    return intersection;
}


COLOR De_que_color(POINT e, POINT d) {
    COLOR color;
    INTERSECTION intersection;
    intersection = First_Intersection(e, d);

    if (intersection.t == INF)
        color = background;
    else {
        SPHERE *obj = (SPHERE *) intersection.object;
        color = obj->color;
    }

    return color;
}

void raytracer() {
    int x_min = viewport.pmin.x;
    int y_min = viewport.pmin.y;
    int x_max = viewport.pmax.x;
    int y_max = viewport.pmax.y;
    POINT w;
    POINT d;
    for (int i = 0; i < H_SIZE; i++) {
        for (int j = 0; j < V_SIZE; j++) {
            w.x = (long double) (i + 1 / 2) * (x_max - x_min) / H_SIZE + x_min;
            w.y = (long double) (j + 1 / 2) * (y_max - y_min) / V_SIZE + y_min;
            w.z = 0;

            long double L = (long double) sqrt(pow(w.x - eye.x, 2) + pow(w.y - eye.y, 2) + pow(w.z - eye.z, 2));
            d.x = (long double) (w.x - eye.x);
            d.y = (long double) (w.y - eye.y);
            d.z = (long double) (w.z - eye.z);
            COLOR color = De_que_color(eye, d);
            plot(i, j, color);
        }
    }
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
    raytracer();
    draw_scene();
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
    // enter GLUT event processing cycle
    glutMainLoop();

    printf("Fin del programa %s...\n\n", argv[0]);


    return 1;
}