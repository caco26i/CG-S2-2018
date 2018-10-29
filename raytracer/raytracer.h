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
} INTERSECTION;

typedef struct {
    float radius;
    POINT center;
} SPHERE;

typedef struct {
	POINT center;
    POINT normal;
    float D;
} PLANE;

typedef struct {
    int n_points;
    POINT** points;
    PLANE* plane;
    int cases;
} POLYGON;


typedef struct {
	INTERSECTION (*fun_ptr)(void*,POINT,POINT);
    POINT (*norm_ptr)(void*,POINT);
    COLOR color;
    float Ka;
    float Kd;
    float Ks;
    float Kn;
    void* object;
} OBJ;

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
OBJ* objects;
VIEWPORT viewport;
LIGHT **lights;

float AmbientIlluminationIntensity;

int window;
int H_SIZE, V_SIZE;


void draw_scene();
void renderScene();
