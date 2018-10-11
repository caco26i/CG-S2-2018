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

void save_file() {
    FILE *fp = fopen("scene.txt", "w");
    if (fp == NULL) return 0;

    for (int i = 0; i < N_LIGHTS; i++){
	    fprintf(fp, "1;%f;%f;%f;%f\n", lights[i].pos.x, lights[i].pos.y, lights[i].pos.z, lights[i].intensity);
    }

	for (i = 0; i < N_SPHERES; i++) {
		fprintf(fp, "2;%f;%f;%f;%f;%f;%f;%f;%f;%f\n", spheres[i].radius, spheres[i].center.x, spheres[i].center.y, spheres[i].center.z, spheres[i].color.R, spheres[i].color.G, spheres[i].color.B, spheres[i].Ka, spheres[i].Kd);

	}
}

int write_truecolor_tga() {
    FILE *fp = fopen("out.tga", "w");
    if (fp == NULL) return 0;

// The image header
    char header[18] = {0}; // char = byte
    header[2] = 2; // truecolor
    header[12] = H_SIZE & 0xFF;
    header[13] = (H_SIZE >> 8) & 0xFF;
    header[14] = V_SIZE & 0xFF;
    header[15] = (V_SIZE >> 8) & 0xFF;
    header[16] = 24; // bits per pixel

    fwrite((const char *) &header, 1, sizeof(header), fp);

// The image data is stored bottom-to-top, left-to-right
    for (int y = V_SIZE - 1; y >= 0; y--)
        for (int x = 0; x < H_SIZE; x++) {
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


    int NX = 2000;
    int NY = 2000;
    int N = 10000000;
    int SCALE = (NX / 8);

    int n,ix,iy;
    double x=0.2,y=0.3,x1=0,y1=0,r=sqrt(3);
    double a0,b0,f1x,f1y;

    for (n=0;n<N_SPHERES;n++) {
        if ((n % (N/10)) == 0)
            fprintf(stderr,".");
        a0 = 3 * (1 + r - x) / (pow(1 + r - x,2.0) + y*y) - (1 + r) / (2 + r);
        b0 = 3 * y / (pow(1 + r - x,2.0) + y*y);
        f1x =  a0 / (a0*a0 + b0*b0);
        f1y = -b0 / (a0*a0 + b0*b0);
        switch (rand()%3) {
            case 0:
                x1 = a0;
                y1 = b0;
                break;
            case 1:
                x1 = -f1x / 2 - f1y * r / 2;
                y1 = f1x * r / 2 - f1y / 2;
                break;
            case 2:
                x1 = -f1x / 2 + f1y * r / 2;
                y1 = -f1x * r / 2 - f1y / 2;
                break;
        }
        if (n < 100)
            continue;
        ix = x * SCALE + NX/2;
        iy = y * SCALE + NY/2;
        x = x1;
        y = y1;
        if (ix < 0 || iy < 0 || ix >= NX || iy >= NY)
            continue;

        float t_x = x*55+H_SIZE/2;
        float t_y =  y*55+V_SIZE/2;
        float dis = sqrt(pow((t_x - (V_SIZE/2)), 2) + pow((t_y - (H_SIZE/2)), 2));
        float dist = sqrt(pow((x), 2) + pow((y), 2));
        spheres[n]->type = T_SPHERE;
        spheres[n]->radius = dist*5+10;
        spheres[n]->center.x = t_x;
        spheres[n]->center.y = t_y;
        spheres[n]->center.z = -dist*20;
        spheres[n]->color.R = (float)((int)(dist*80)%255) / 255;
        spheres[n]->color.G = (float)((int)(dist*20)%255) / 255;
        spheres[n]->color.B = 0.5;
        spheres[n]->Ka = 1;
        spheres[n]->Kd = 1;

		printf("%f\n", (float)spheres[n]->center.z);
    }

//
//    for (int i = 0; i < N_SPHERES; ++i)
//    {
//        spheres[i]->type = T_SPHERE;
//        spheres[i]->radius = (rand()%20)+10;
//        spheres[i]->center.x = (rand()%1000);
//        spheres[i]->center.y = (rand()%1000);
//        spheres[i]->center.z = (rand()%100)+100;
//        spheres[i]->color.R = (rand()%100)/100.0;
//        spheres[i]->color.G = (rand()%100)/100.0;
//        spheres[i]->color.B = (rand()%100)/100.0;
//        spheres[i]->Ka = 1;
//        spheres[i]->Kd = 1;
//    }
}

void initLights(){
    lights = malloc(sizeof(LIGHTSOURCE)* N_LIGHTS);
    for (int i = 0; i < N_LIGHTS; ++i)
    {
        lights[i].pos.x = H_SIZE/2;
        lights[i].pos.y = V_SIZE/2;
        lights[i].pos.z = -150;
        lights[i].intensity = 1;
    }
}

void init() {
    srand (time(NULL));

    AmbientIlluminationIntensity = 0.2;

    eye.x = H_SIZE/ 2;
    eye.y =  V_SIZE/ 2;
    eye.z = -1000;

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
    initLights();
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
    //float a = pow((d.x - e.x), 2.0) + pow((d.y - e.y), 2.0) + pow((d.z - e.z), 2.0);
    float a = 1;
    float b = 2.0 * ((d.x )*( e.x - sphere->center.x) + (d.y )*( e.y - sphere->center.y) +
                     (d.z )*( e.z - sphere->center.z));
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


COLOR De_que_color(POINT e, POINT d, int i, int j) {
    COLOR color;
    INTERSECTION intersection;
    intersection = First_Intersection(e, d);

    if (intersection.t == INF){
        color.R = (float)j/1000.0;
        color.G = (float)j/1000.0;
        color.B = 1.0;
    }

    else {
        SPHERE *obj = (SPHERE *) intersection.object;
        color = obj->color;

        POINT intersectionPoint;
        float intensity = 0;

        intersectionPoint.x = e.x + intersection.t * d.x;
        intersectionPoint.y = e.y + intersection.t * d.y;
        intersectionPoint.z = e.z + intersection.t * d.z;
        float n;
        POINT N;
        N.x = intersectionPoint.x - obj->center.x;
        N.y = intersectionPoint.y - obj->center.y;
        N.z = intersectionPoint.z - obj->center.z;
        n = sqrt(pow(N.x, 2) + pow(N.y, 2) + pow(N.z, 2));
        N.x /=n;
        N.y /=n;
        N.z /=n;

        POINT L;

        for (int i = 0; i < N_LIGHTS; ++i)
        {

            L.x = lights[i].pos.x - intersectionPoint.x;
            L.y = lights[i].pos.y - intersectionPoint.y;
            L.z = lights[i].pos.z - intersectionPoint.z;
            n = sqrt(pow(L.x, 2) + pow(L.y, 2) + pow(L.z, 2));
            L.x /=n;
            L.y /=n;
            L.z /=n;


            intensity +=((L.x * N.x + L.y * N.y + L.z * N.z) * obj->Kd * lights[i].intensity);
        }
        intensity += obj->Ka * AmbientIlluminationIntensity;
        if(intensity<0.0)intensity=0.0;
        if(intensity>1.0)intensity=1.0;

        color.R *=intensity;
        color.G *=intensity;
        color.B *=intensity;
    }

    return color;
}

void raytracer() {
    int x_min = viewport.pmin.x;
    int y_min = viewport.pmin.y;
    int x_max = viewport.pmax.x;
    int y_max = viewport.pmax.y;
    POINT d;
    POINT w;
    for (int i = 0; i < H_SIZE; i++) {
        for (int j = 0; j < V_SIZE; j++) {
            w.x = i + 0.5;//((i * (x_max - x_min)) / H_SIZE) + x_min;
            w.y = j + 0.5;//((j * (y_max - y_min)) / V_SIZE) + y_min;
            w.z = 0;

            float L = sqrt(pow(w.x - eye.x, 2) + pow(w.y - eye.y, 2) + pow(w.z - eye.z, 2));
            d.x = (w.x-eye.x)/L;
            d.y = (w.y-eye.y)/L;
            d.z = (w.z-eye.z)/L;

            COLOR color = De_que_color(eye, d, i, j);
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
