/*
 * Instituto Tecnologico de Costa Rica
 * Escuela de Ingenieria en Computacion
 * Computer Graphics
 *
 * Programa: Proyecto 1
 * Archivo:  proyecto_1.h
 */

typedef struct {
    float x;
    float y;
    float z;
} POINT;

typedef struct {
    float R;
    float G;
    float B;
} COLOR;


typedef struct {
    float t;
    void* object;
    bool set;
} INTERSECTION;

typedef struct {
    int type;
    INTERSECTION (*fun_ptr)(void*,POINT,POINT); 
    float radius;
    POINT center;
    COLOR color;
    float Ka;
    float Kd;
    float Ks;
    float Kn;
} SPHERE;


typedef struct {
    POINT pos;
    float intensity;
} LIGHT;

typedef struct {
    POINT pmin;
    POINT pmax;
} VIEWPORT;

POINT eye;
COLOR **buffer;
COLOR color;
COLOR background;
void**objects;
VIEWPORT viewport;
LIGHT **lights;

float AmbientIlluminationIntensity;

int window;
int H_SIZE, V_SIZE;


void draw_scene();
void renderScene();
