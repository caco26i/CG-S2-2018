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
    POINT center;
    COLOR color;
    float Ka;
    float Kd;
    float Ks;
    float Kn;
} OBJ;

typedef struct {
    int type;
    float radius;
    OBJ object;
} SPHERE;

typedef struct {
    int type;
    POINT normal;
    OBJ object;
    float D;
} PLANE;

typedef struct {
    int n_points;
    POINT** points;
    PLANE* plane;
    int cases;
} POLYGON;

typedef struct {
    float t;
    OBJ object;
    COLOR color;
    bool set;
} INTERSECTION;

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
SPHERE **spheres;
POLYGON **polygons;
VIEWPORT viewport;
LIGHT **lights;

float AmbientIlluminationIntensity;

int window;
int H_SIZE, V_SIZE;


void draw_scene();
void renderScene();
