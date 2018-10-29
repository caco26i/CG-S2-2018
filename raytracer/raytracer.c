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
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "raytracer.h"

#define N_RAYS 4
#define INF 10000000
#define SHADOW_K 1
#define PROGRESS 0
#define DEBUG 0

int H_SIZE;
int V_SIZE;
int N_POLYGONS;
int N_SPHERES;
int N_LIGHTS;
int ANTIALIASING;
int SHADOWS;
int N_OBJECTS;
int N_PLANES;

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

void set_color(float r, float g, float b) {
    color.R = r;
    color.G = g;
    color.B = b;
}

float myPow(float num, int exp){
	float res = num;
	for (int i = 1; i < exp; ++i)
	{
		res *= num;
	}
	return res;
}


POINT dot(POINT a, POINT b) {
    POINT dot;
    dot.x = (a.y*b.z) - (a.z*b.y);
    printf("DOT.X %f\n", (a.y*b.z));
    dot.y = -((a.x*b.z) - (a.z*b.x));
    dot.z = (a.x*b.y) - (a.y*b.x);
    return dot;
}

POINT T(POINT P, int Dx, int Dy) {
    P.x += Dx;
    P.y += Dy;
    return P;
}

bool esta_punto_poly(POLYGON* poly, POINT punto){
    float X;
    float Y;
    float testx;
    float testy;
    POINT points[poly->n_points];

//    testx = punto.x;
//    testy = punto.y;

    for (int i = 0; i < poly->n_points; ++i){
        if (poly->cases == 0){
            points[i].x = poly->points[i]->y;
            points[i].y = poly->points[i]->z;
            testx = punto.y;
            testy = punto.z;
        } else if (poly->cases == 1) {
            points[i].x = poly->points[i]->x;
            points[i].y = poly->points[i]->z;
            testx = punto.x;
            testy = punto.z;
        }
        else {
            points[i].x = poly->points[i]->x;
            points[i].y = poly->points[i]->y;
            testx = punto.x;
            testy = punto.y;
        }
    }

    int i, j, c = 0;
    for (i = 0, j = poly->n_points-1; i < poly->n_points; j = i++) {
        if ( ((points[i].y>testy) != (points[j].y>testy)) &&
             (testx < (points[j].x-points[i].x) * (testy-points[i].y) / (points[j].y-points[i].y) + points[i].x) ) {
            c = !c;
        }
    }
//    if(c) printf("TEST X %f Y %f Z %f\n", punto.x, punto.y, punto.z);
//    printf("CASES %d\n", poly->cases);
    return c;

}


INTERSECTION IntersectionPoly(void* obj, POINT e, POINT d){
	if(DEBUG)printf("IntersectionPoly\n");
	POLYGON * polygon = (POLYGON *)((OBJ*)obj)->object;

    INTERSECTION intersection;
    PLANE* plane = polygon->plane;

    float denominador = (plane->normal.x * d.x + plane->normal.y * d.y + plane->normal.z * d.z);

    if (denominador == 0){
        intersection.t = INF;
    }else{
        intersection.t =
                -((plane->normal.x * e.x + plane->normal.y * e.y + plane->normal.z * e.z) + plane->D )
                / denominador;

        POINT intersection_point;
        intersection_point.x = e.x + intersection.t*d.x;
        intersection_point.y = e.y + intersection.t*d.y;
        intersection_point.z = e.z + intersection.t*d.z;

        if(intersection.t < 0.0001)intersection.t = INF;
        else if(polygon->n_points > 2 && !esta_punto_poly(polygon, intersection_point)) intersection.t = INF;
//    else {printf("PUNTO X %f Y %f Z %f\n", intersection_point.x, intersection_point.y, intersection_point.z);}
    }

    intersection.object = obj;
    return intersection;
}

INTERSECTION IntersectionSphere(void* obj, POINT e, POINT d){

	if(DEBUG)printf("IntersectionSphere\n");
	SPHERE* sphere = (SPHERE*)((OBJ*)obj)->object;
	//printf("%f \n", sphere->radius);
    INTERSECTION inter;
    //float a = pow((d.x - e.x), 2.0) + pow((d.y - e.y), 2.0) + pow((d.z - e.z), 2.0);
    float a = 1;
    float b = 2.0 * ((d.x )*( e.x - sphere->center.x) + (d.y )*( e.y - sphere->center.y) +
                     (d.z )*( e.z - sphere->center.z));
    float g = myPow((e.x - sphere->center.x), 2.0) + myPow((e.y - sphere->center.y), 2.0) +
              myPow((e.z - sphere->center.z), 2.0) - myPow(sphere->radius, 2.0);
    float delta = myPow(b, 2.0) - 4.0 * g * a;


    if(delta < -0.001){
        inter.t = INF;
    }else if(delta < 0.001){
        inter.t = -b/(2.0*a);

    }else{
        float t1,t2;
        t1 = (-b + sqrt(delta))/(2.0*a);
        t2 = (-b - sqrt(delta))/(2.0*a);
        if(t1 < -0.001 && t2 < -0.001){
            inter.t = INF;
        }else if(t1 < -0.001){
            inter.t = t2;
        }else if(t2 < -0.001){
            inter.t = t1;
        }else{
            inter.t = (t1<t2)?t1:t2;
        }
    }
    if(inter.t < -0.001)inter.t = INF;
    inter.object = obj;
    return inter;
}

INTERSECTION IntersectionPlane(void* obj, POINT e, POINT d){
	if(DEBUG)printf("IntersectionSphere\n");
	PLANE* plane = (PLANE*)((OBJ*)obj)->object;
    INTERSECTION inter;
    long double D = -(plane->normal.x * plane->center.x + plane->normal.y * plane->center.y + plane->normal.z * plane->center.z);

    inter.t =
            ((plane->normal.x * e.x + plane->normal.y * e.y + plane->normal.z * e.z) + D )
            / (plane->normal.x * d.x + plane->normal.y * d.y + plane->normal.z * d.z);

    POINT intersection_point;
    intersection_point.x = e.x + inter.t*d.x;
    intersection_point.y = e.y + inter.t*d.y;
    intersection_point.z = e.z + inter.t*d.z;

    if(inter.t < -0.001)inter.t = INF;
    inter.object = obj;
    return inter;
}

POINT GetNormalSphere(void* obj, POINT i){
	SPHERE* sphere = (SPHERE*)((OBJ*)obj)->object;
	POINT N;
	N.x = i.x - sphere->center.x;
    N.y = i.y - sphere->center.y;
    N.z = i.z - sphere->center.z;
    float n = sqrt(myPow(N.x, 2) + myPow(N.y, 2) + myPow(N.z, 2));
    N.x /=n;
    N.y /=n;
    N.z /=n;
    return N;
}

POINT GetNormalPoly(void* obj, POINT i){

	POLYGON * polygon = (POLYGON *)((OBJ*)obj)->object;
	return polygon->plane->normal;
}

POINT GetNormalPlane(void* obj, POINT i){

	PLANE * plane = (PLANE *)((OBJ*)obj)->object;
	return plane->normal;
}

void loadScene(){
    scanf("H_SIZE %d\n",&H_SIZE);
    scanf("V_SIZE %d\n",&V_SIZE);
    scanf("ANTIALIASING %d\n",&ANTIALIASING);
    scanf("SHADOWS %d\n",&SHADOWS);
    scanf("N_LIGHTS %d\n",&N_LIGHTS);
    scanf("N_SPHERES %d\n",&N_SPHERES);
    scanf("N_PLANES %d\n",&N_PLANES);
    scanf("N_POLYGONS %d\n",&N_POLYGONS);

    printf("H_SIZE %d\n",H_SIZE);
    printf("V_SIZE %d\n",V_SIZE);
    printf("ANTIALIASING %d\n",ANTIALIASING);
    printf("SHADOWS %d\n",SHADOWS);
    printf("N_LIGHTS %d\n",N_LIGHTS);
    printf("N_SPHERES %d\n",N_SPHERES);
    printf("N_PLANES %d\n",N_PLANES);
    printf("N_POLYGONS %d\n",N_POLYGONS);


    N_OBJECTS = N_POLYGONS + N_SPHERES + N_PLANES;
    printf("N_OBJECTS %d\n",N_OBJECTS);
    int i = 0;
    int j;
    objects = malloc(sizeof(OBJ)*N_OBJECTS);


    j = i;
    for (i; i < N_SPHERES + j; ++i)
    {
        objects[i].object = (void*)malloc(sizeof(SPHERE));
        SPHERE* auxSphere = (SPHERE*)objects[i].object;
        scanf("radius %f\n",&auxSphere->radius);
        scanf("x %f y %f z %f\n",&auxSphere->center.x,&auxSphere->center.y,&auxSphere->center.z);
        scanf("R %f G %f B %f\n",&objects[i].color.R,&objects[i].color.G,&objects[i].color.B);
        //printf("R %f G %f B %f\n",objects[i].color.R,objects[i].color.G,objects[i].color.B);
        scanf("Ka %f Kd %f Ks %f Kn %f\n",&objects[i].Ka,&objects[i].Kd,&objects[i].Ks,&objects[i].Kn);
        objects[i].fun_ptr = &IntersectionSphere;
        objects[i].norm_ptr = &GetNormalSphere;
        //objects[i].type = 0;
    }
    j = i;
    for (i; i < N_POLYGONS + j; ++i) {

        int n_points;
        scanf("N_POINTS %d\n",&n_points);
        printf("N_POINTS %d\n",n_points);

        objects[i].object = (void*)malloc(sizeof(POLYGON));
        POLYGON* auxPolygon = (POLYGON*)objects[i].object;
        auxPolygon->points = malloc(sizeof(POINT)*n_points);
        auxPolygon->n_points = n_points;

        for (int k = 0; k < n_points; ++k)
        {
            auxPolygon->points[k] = malloc(sizeof(POINT));
            scanf("x %f y %f z %f\n", &auxPolygon->points[k]->x, &auxPolygon->points[k]->y, &auxPolygon->points[k]->z);
            printf("x %f y %f z %f\n", auxPolygon->points[k]->x, auxPolygon->points[k]->y, auxPolygon->points[k]->z);

        }

        auxPolygon->plane = malloc(sizeof(PLANE));

        //Se crean los vectores de direccion
        if(n_points >= 3){

            POINT A;
            A.x = auxPolygon->points[1]->x - auxPolygon->points[0]->x;
            A.y = auxPolygon->points[1]->y - auxPolygon->points[0]->y;
            A.z = auxPolygon->points[1]->z - auxPolygon->points[0]->z;
            POINT B;
            B.x = auxPolygon->points[2]->x - auxPolygon->points[0]->x;
            B.y = auxPolygon->points[2]->y - auxPolygon->points[0]->y;
            B.z = auxPolygon->points[2]->z - auxPolygon->points[0]->z;

            auxPolygon->plane->normal = dot(A, B);
        }
        else {
            auxPolygon->plane->normal.x = auxPolygon->points[1]->x;
            auxPolygon->plane->normal.y = auxPolygon->points[1]->y;
            auxPolygon->plane->normal.z = auxPolygon->points[1]->z;
        }

        auxPolygon->plane->D =
                -((auxPolygon->plane->normal.x * auxPolygon->points[0]->x)
                  + auxPolygon->plane->normal.y * auxPolygon->points[0]->y
                  +(auxPolygon->plane->normal.z * auxPolygon->points[0]->z));

        objects[i].fun_ptr = &IntersectionPoly;
        objects[i].norm_ptr = &GetNormalPoly;
        auxPolygon->plane->center.x = auxPolygon->points[0]->x;
        auxPolygon->plane->center.y = auxPolygon->points[0]->y;
        auxPolygon->plane->center.z = auxPolygon->points[0]->z;


        if (fabs(auxPolygon->plane->normal.x) > fabs (auxPolygon->plane->normal.y) && fabs(auxPolygon->plane->normal.x) > fabs(auxPolygon->plane->normal.z)){
            auxPolygon->cases = 0;
        } else if (fabs (auxPolygon->plane->normal.y) > fabs (auxPolygon->plane->normal.x) && fabs (auxPolygon->plane->normal.y) > fabs (auxPolygon->plane->normal.z) ){
            auxPolygon->cases = 1;
        } else if (fabs (auxPolygon->plane->normal.z) > fabs (auxPolygon->plane->normal.x) && fabs (auxPolygon->plane->normal.z) > fabs (auxPolygon->plane->normal.y) ){
            auxPolygon->cases = 2;
        }

        scanf("R %f G %f B %f\n",&objects[i].color.R,&objects[i].color.G,&objects[i].color.B);
        scanf("Ka %f Kd %f Ks %f Kn %f\n",&objects[i].Ka,&objects[i].Kd,&objects[i].Ks,&objects[i].Kn);
    }
    j = i;
    for (i; i < N_PLANES + j; ++i) {

        objects[i].object = (void*)malloc(sizeof(PLANE));
        PLANE* auxPlane = (PLANE*)objects[i].object;
        POINT * points = malloc(sizeof(POINT)*3);

        for (int k = 0; k < 3; ++k)
        {
            scanf("x %f y %f z %f\n", &points[k].x, &points[k].y, &points[k].z);
        }

        //Se crean los vectores de direccion
        POINT A;
        A.x = points[1].x - points[0].x;
        A.y = points[1].y - points[0].y;
        A.z = points[1].z - points[0].z;
        POINT B;
        B.x = points[2].x - points[0].x;
        B.y = points[2].y - points[0].y;
        B.z = points[2].z - points[0].z;

        objects[i].fun_ptr = &IntersectionPlane;
        objects[i].norm_ptr = &GetNormalPlane;
        auxPlane->normal = dot(A, B);
        auxPlane->center.x = points[0].x;
        auxPlane->center.y = points[0].y;
        auxPlane->center.z = points[0].z;

        scanf("R %f G %f B %f\n",&objects[i].color.R,&objects[i].color.G,&objects[i].color.B);
        scanf("Ka %f Kd %f Ks %f Kn %f\n",&objects[i].Ka,&objects[i].Kd,&objects[i].Ks,&objects[i].Kn);
    }
    lights = malloc(sizeof(LIGHT*)*N_LIGHTS);

    for (int i = 0; i < N_LIGHTS; ++i)
    {
        lights[i] = malloc(sizeof(LIGHT));
        scanf("x %f y %f z %f\n",&lights[i]->pos.x,&lights[i]->pos.y,&lights[i]->pos.z);
        scanf("intensity %f\n",&lights[i]->intensity);

    }
}

void init() {
    srand (time(NULL));

    AmbientIlluminationIntensity = 0.2;

    eye.x = H_SIZE/2;
    eye.y =  V_SIZE/2;
    eye.z = -1000;

    color.R = 0;
    color.G = 0;
    color.B = 0;

    background.R = 0.3;
    background.G = 0;
    background.B = 0.1;

    viewport.pmin.x = 0;
    viewport.pmin.y = 0;
    viewport.pmin.z = 0;
    viewport.pmax.x = H_SIZE;
    viewport.pmax.y = V_SIZE;
    viewport.pmax.z = 0;
}

INTERSECTION First_Intersection(POINT e, POINT d) {
    if(DEBUG)printf("First_Intersection\n");
    float tmin;
    INTERSECTION intersection;
    intersection.t = INF;
    INTERSECTION auxIntersection;

    for (int i = 0; i < N_OBJECTS; ++i)
    {
        if(DEBUG)printf("------First_Intersection %d\n", i);

        auxIntersection = (objects[i].fun_ptr)((void*)&objects[i], e, d);

        if(auxIntersection.t > (float)SHADOW_K && auxIntersection.t < intersection.t){
            intersection.t = auxIntersection.t;
            intersection.object = auxIntersection.object;
        }

    }
    return intersection;
}

COLOR De_que_color(POINT e, POINT d) {
	INTERSECTION intersection;
	//if(DEBUG)printf("De_que_color\n");
    intersection = First_Intersection(e, d);
    COLOR color;

    if (intersection.t == INF){
        color = background;
    }

    else {


    	//printf("pega en algo\n");
        OBJ* obj = (OBJ*)intersection.object;
		color = ((OBJ*)intersection.object)->color;
		//printf("interseccion R %f G %f B %f\n",color.R,color.G,color.B);
        POINT intersectionPoint;
        float intensity = 0;
        float E = 0;

        intersectionPoint.x = e.x + intersection.t * d.x;
        intersectionPoint.y = e.y + intersection.t * d.y;
        intersectionPoint.z = e.z + intersection.t * d.z;
        float n;
        POINT N = (obj->norm_ptr)((void*)obj, intersectionPoint);
        POINT L;

        for (int i = 0; i < N_LIGHTS; ++i)
        {


    		float Fatt,cosNL,cosVR;

            L.x = lights[i]->pos.x - intersectionPoint.x;
            L.y = lights[i]->pos.y - intersectionPoint.y;
            L.z = lights[i]->pos.z - intersectionPoint.z;
            n = sqrt(myPow(L.x, 2) + myPow(L.y, 2) + myPow(L.z, 2));
            L.x /=n;
            L.y /=n;
            L.z /=n;


            bool lid = true;
            if(SHADOWS){
            	INTERSECTION intersectionLight;
            	intersectionLight = First_Intersection(intersectionPoint, L);

	            if (intersectionLight.t != INF && intersectionLight.t < n){
	            	lid = false;
			    }
            }



            if(lid){


            	POINT R;
            	POINT V;

            	Fatt = fmin(1.0,1.0 / (myPow(n/(lights[i]->intensity * 1000.0),2)));

            	cosNL = (L.x * N.x + L.y * N.y + L.z * N.z);

            	R.x = 2 * N.x * cosNL - L.x;
            	R.y = 2 * N.y * cosNL - L.y;
            	R.z = 2 * N.z * cosNL - L.z;

            	V.x = -d.x;
            	V.y = -d.y;
            	V.z = -d.z;

            	cosVR = (V.x * R.x + V.y * R.y + V.z * R.z);


            	if(cosNL > 0){
            		if(cosVR > 0)E += (myPow(cosVR,obj->Kn) * obj->Ks * lights[i]->intensity * Fatt);
            		intensity += (cosNL * obj->Kd * lights[i]->intensity * Fatt);
            	}


            }




        }

        intensity += obj->Ka * AmbientIlluminationIntensity;

        if(intensity>1.0)intensity=1.0;
		if(E>1.0)E=1.0;

        color.R *=intensity;
        color.G *=intensity;
        color.B *=intensity;

        color.R = color.R + E * (1 - color.R);
        color.G = color.G + E * (1 - color.G);
        color.B = color.B + E * (1 - color.B);


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
    if(DEBUG)printf("raytracer\n");
    for (int i = 0; i < H_SIZE; i++) {
    	if(DEBUG)printf("i = %d\n",i);
    	if(PROGRESS)printf("%f\n", ((float)i / (float)H_SIZE)*100);
        for (int j = 0; j < V_SIZE; j++) {
        	//if(DEBUG)printf("j = %d\n",j);

        	COLOR color;
        	color.R = 0;
        	color.G = 0;
        	color.B = 0;

        	if(!ANTIALIASING){
        		w.x = i + 0.5;//((i * (x_max - x_min)) / H_SIZE) + x_min;
	            w.y = j + 0.5;//((j * (y_max - y_min)) / V_SIZE) + y_min;
	            w.z = 0;

	            float L = sqrt(myPow(w.x - eye.x, 2) + myPow(w.y - eye.y, 2) + myPow(w.z - eye.z, 2));
	            d.x = (w.x-eye.x)/L;
	            d.y = (w.y-eye.y)/L;
	            d.z = (w.z-eye.z)/L;

	            color = De_que_color(eye, d);
        	}else{
        		for (int k = 0; k < N_RAYS / 2; ++k)
	        	{
	        		for (int l = 0; l < N_RAYS / 2; ++l)
	        		{
	        			w.x = i + (1/(N_RAYS / 2 - 1))*k;//((i * (x_max - x_min)) / H_SIZE) + x_min;
			            w.y = j + (1/(N_RAYS / 2 - 1))*l;//((j * (y_max - y_min)) / V_SIZE) + y_min;
			            w.z = 0;

			            float L = sqrt(myPow(w.x - eye.x, 2) + myPow(w.y - eye.y, 2) + myPow(w.z - eye.z, 2));
			            d.x = (w.x-eye.x)/L;
			            d.y = (w.y-eye.y)/L;
			            d.z = (w.z-eye.z)/L;

			            COLOR auxColor = De_que_color(eye, d);
			            color.R += auxColor.R;
			        	color.G += auxColor.G;
			        	color.B += auxColor.B;
	        		}

	        	}

	            color.R /= N_RAYS;
	        	color.G /= N_RAYS;
	        	color.B /= N_RAYS;
        	}

            plot(i, j, color);
        }
    }
    printf("Listo\n");
}

int main(int argc, char **argv) {
    int i, j, length;

    loadScene();
    init();
    buffer = (COLOR **) malloc(H_SIZE * sizeof(COLOR *));
    for (i = 0; i < H_SIZE; i++) {
        buffer[i] = (COLOR *) malloc(V_SIZE * sizeof(COLOR));
    }

    set_color(1, 1, 1);


    raytracer();
    write_truecolor_tga();

    system("convert out.tga out.png");
    printf("Fin del programa %s...\n\n", argv[0]);


    return 1;
}
