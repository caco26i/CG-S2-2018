/*
 * Instituto Tecnologico de Costa Rica
 * Escuela de Ingenieria en Computacion
 * Computer Graphics
 *
 * Programa: Proyecto 1
 * Archivo:  proyecto_1.h
 */
typedef struct {
    double x;
    double y;
    double z;
} POINT;

typedef struct {
    double R;
    double G;
    double B;
} COLOR;


typedef struct {
    double t;
    double m;
    POINT collision;
    void* object;
} INTERSECTION;

typedef struct { // type 1
    double radius;
    POINT center;
} SPHERE;

typedef struct { // type 2
	POINT center;
    POINT normal;
    double D;
} PLANE;

typedef struct { // type 3
    int n_points;
    POINT** points;
    PLANE* plane;
    int cases;
} POLYGON;

typedef struct { // type 4
    double radius;
    double d1;
    double d2;
    POINT center;
    POINT axis;
} CYLINDER;

typedef struct { // type 5
    double radius1;
    double radius2;
    PLANE* plane;
} DISC;

typedef struct { // type 6
    double radius;
    POINT center1;
    POINT center2;
    POINT normal;
} OVAL;


typedef struct {
	INTERSECTION (*fun_ptr)(void*,POINT,POINT);
    POINT (*norm_ptr)(INTERSECTION);
    COLOR color;
    double Ka;
    double Kd;
    double Ks;
    double Kn;
    void* object;
    double O1;
    double O2;
    double O3;
} OBJ;

typedef struct {
    POINT pos;
    double intensity;
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

double AmbientIlluminationIntensity;

int window;
int H_SIZE, V_SIZE;


void draw_scene();
void renderScene();
