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
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "raytracer.h"

int write_truecolor_tga() {
    FILE *fp = fopen("out.tga", "w");
    if (fp == NULL) return 0;

// The image header
    char header[18] = {0}; // char = byte
    header[2] = 2; // truecolor
    header[12] = V_SIZE & 0xFF;
    header[13] = (V_SIZE >> 8) & 0xFF;
    header[14] = H_SIZE & 0xFF;
    header[15] = (H_SIZE >> 8) & 0xFF;
    header[16] = 24; // bits per pixel

    fwrite((const char *) &header, 1, sizeof(header), fp);

// The image data is stored bottom-to-top, left-to-right
    for (int y = H_SIZE - 1; y >= 0; y--)
        for (int x = 0; x < V_SIZE; x++) {
            fputc((int) (buffer[x][y].B * 255), fp);
            fputc((int) (buffer[x][y].G * 255), fp);
            fputc((int) (buffer[x][y].R * 255), fp);
        }

// The file footer
    static const char footer[26] =
            "\0\0\0\0" // no extension area
            "\0\0\0\0" // no developer directory
            "TRUEVISION-XFILE" // yep, this is a TGA file
            ".";
    fwrite((const char *) &footer, 1, sizeof(footer), fp);

    fclose(fp);
    return 1;
}

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
        spheres[i]->radius = (rand()%100)+10;
        spheres[i]->center.x = (rand()%1000);
        spheres[i]->center.y = (rand()%1000);
        spheres[i]->center.z = (rand()%100)+100;
        spheres[i]->color.R = (rand()%100)/100.0;
        spheres[i]->color.G = (rand()%100)/100.0;
        spheres[i]->color.B = (rand()%100)/100.0;

    }
}

void init() {
    srand (time(NULL));

    eye.x = H_SIZE / 2;
    eye.y = V_SIZE / 2;
    eye.z = -100;

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
    viewport.pmax.x = H_SIZE;
    viewport.pmax.y = V_SIZE;
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
    float a = pow((d.x - e.x), 2.0) + pow((d.y - e.y), 2.0) + pow((d.z - e.z), 2.0);
    float b = 2.0 * ((d.x - e.x) * (e.x - sphere->center.x) + (d.y - e.y) * (e.y - sphere->center.y) +
                     (d.z - e.z) * (e.z - sphere->center.z));
    float g = pow((e.x - sphere->center.x), 2.0) + pow((e.y - sphere->center.y), 2.0) +
              pow((e.z - sphere->center.z), 2.0) - pow(sphere->radius, 2.0);
    float delta = pow(b, 2.0) - 4.0 * g * a;


    if (delta > 0.0) {

    } else {
        /*
        printf("pego\n");
        printf("a %f\n", a);
        printf("b %f\n", b);
        printf("g %f\n", g);
        printf("pow(b,2) %f\n", pow(b,2));
        printf("4.0 * g * a %f\n", 4.0 * g * a);*/
    }

    if(delta < -0.001){
        intersection.t = INF;
    }else if(delta < 0.001){
        intersection.t = -b/(2.0*a);

    }else{
        float t1,t2;
        t1 = (-b + sqrt(delta))/(2.0*a);
        t2 = (-b - sqrt(delta))/(2.0*a);
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
    POINT d;
    for (int i = 0; i < H_SIZE; i++) {
        for (int j = 0; j < V_SIZE; j++) {
            d.x = i;//((i * (x_max - x_min)) / H_SIZE) + x_min;
            d.y = j;//((j * (y_max - y_min)) / V_SIZE) + y_min;
            d.z = 0;

            //float L = sqrt(pow(w.x - eye.x, 2) + pow(w.y - eye.y, 2) + pow(w.z - eye.z, 2));

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
    write_truecolor_tga();
}

int main(int argc, char **argv) {
    int i, j, length;

//    printf("Resolución: %s\n", argv[1]);
//    length = strlen(argv[1]);
//    for (i = 0; i < length; i++)
//        if (!isdigit(argv[1][i])) {
//            printf("Ingrese una resulución válida\n");
//            return 1;
//        }

    // Set Res
    H_SIZE = 1008;
    V_SIZE = 567;
//    sscanf(argv[1], "%d", &H_SIZE);
//    sscanf(argv[1], "%d", &V_SIZE);

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