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
#define DEBUG 0



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

double myPow(double num, int exp){
	double res = num;
	for (int i = 1; i < exp; ++i)
	{
		res *= num;
	}
	return res;
}


POINT dot(POINT a, POINT b) {
    POINT dot;
    dot.x = (a.y*b.z) - (a.z*b.y);
    //printf("DOT.X %lf\n", (a.y*b.z));
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
    double X;
    double Y;
    double testx;
    double testy;
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
//    if(c) printf("TEST X %lf Y %lf Z %lf\n", punto.x, punto.y, punto.z);
//    printf("CASES %d\n", poly->cases);
    return c;

}


INTERSECTION IntersectionPoly(void* obj, POINT e, POINT d){
	if(DEBUG)printf("IntersectionPoly\n");
	POLYGON * polygon = (POLYGON *)((OBJ*)obj)->object;

    INTERSECTION inter;
    PLANE* plane = polygon->plane;

    double denominador = (plane->normal.x * d.x + plane->normal.y * d.y + plane->normal.z * d.z);

    if (denominador == 0){
        inter.t = INF;
    }else{
        inter.t =
                -((plane->normal.x * e.x + plane->normal.y * e.y + plane->normal.z * e.z) + plane->D )
                / denominador;

        POINT intersection_point;
        intersection_point.x = e.x + inter.t*d.x;
        intersection_point.y = e.y + inter.t*d.y;
        intersection_point.z = e.z + inter.t*d.z;

        if(inter.t < 0.0001)inter.t = INF;
        else if(polygon->n_points > 2 && !esta_punto_poly(polygon, intersection_point)) inter.t = INF;
//    else {printf("PUNTO X %lf Y %lf Z %lf\n", intersection_point.x, intersection_point.y, intersection_point.z);}
    }

    if(inter.t != INF){
        //printf("Colision con Cilindro\n");
        inter.collision.x = e.x + inter.t * d.x;
        inter.collision.y = e.y + inter.t * d.y;
        inter.collision.z = e.z + inter.t * d.z;
        inter.object = obj;
    }
    return inter;
}

INTERSECTION IntersectionSphere(void* obj, POINT e, POINT d){

	if(DEBUG)printf("IntersectionSphere\n");
	SPHERE* sphere = (SPHERE*)((OBJ*)obj)->object;
	//printf("%lf \n", sphere->radius);
    INTERSECTION inter;
    //double a = pow((d.x - e.x), 2.0) + pow((d.y - e.y), 2.0) + pow((d.z - e.z), 2.0);

    POINT X;
    X.x = e.x - sphere->center.x;
    X.y = e.y - sphere->center.y;
    X.z = e.z - sphere->center.z;

    double a = 1;
    double b = 2.0 * ((d.x )*( X.x) + (d.y )*( X.y) +
                     (d.z )*( X.z));
    double g = myPow((X.x), 2.0) + myPow((X.y), 2.0) +
              myPow((X.z), 2.0) - myPow(sphere->radius, 2.0);
    double delta = myPow(b, 2.0) - 4.0 * g * a;


    if(delta < -0.001){
        inter.t = INF;
    }else if(delta < 0.001){
        inter.t = -b/(2.0*a);

    }else{
        double t1,t2;
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
    if(inter.t < 0.001)inter.t = INF;
    if(inter.t != INF){
        //printf("Colision con Cilindro\n");
        inter.collision.x = e.x + inter.t * d.x;
        inter.collision.y = e.y + inter.t * d.y;
        inter.collision.z = e.z + inter.t * d.z;
        inter.object = obj;
    }
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

    
    if(inter.t < -0.001)inter.t = INF;
    if(inter.t != INF){
        //printf("Colision con Cilindro\n");
        inter.collision.x = e.x + inter.t * d.x;
        inter.collision.y = e.y + inter.t * d.y;
        inter.collision.z = e.z + inter.t * d.z;
        inter.object = obj;
    }
    
    return inter;
}

INTERSECTION IntersectionCylinder(void* obj, POINT e, POINT d){

    if(DEBUG)printf("IntersectionCylinder\n");
    CYLINDER* cylinder = (CYLINDER*)((OBJ*)obj)->object;
    
    INTERSECTION inter;

    POINT X;
    X.x = e.x - cylinder->center.x;
    X.y = e.y - cylinder->center.y;
    X.z = e.z - cylinder->center.z;

    inter.t == INF;
    double t1,t2;
    double a = myPow(d.x,2) + myPow(d.y,2) + myPow(d.z,2) - myPow((d.x*cylinder->axis.x + d.y*cylinder->axis.y + d.z*cylinder->axis.z),2);

    double b = 2.0 * (
        ((d.x )*(X.x) + (d.y )*(X.y) + (d.z)*(X.z)) - 
        (d.x*cylinder->axis.x + d.y*cylinder->axis.y + d.z*cylinder->axis.z) * 
        ((cylinder->axis.x )*(X.x) + (cylinder->axis.y)*(X.y) + (cylinder->axis.z)*(X.z)) );

    double g = 
            myPow((X.x), 2.0) + myPow((X.y), 2.0) + myPow((X.z), 2.0) - 
            myPow(((cylinder->axis.x )*(X.x) + (cylinder->axis.y)*(X.y) + (cylinder->axis.z)*(X.z)), 2.0) - 
            myPow(cylinder->radius, 2.0);

    double delta = myPow(b, 2.0) - 4.0 * g * a;

    //printf("%lf\n",delta );
    if(delta < -0.001){
        //printf("no hay Colision con Cilindro\n");
        t1 = INF;
        t2 = INF;
    }else if(delta < 0.001){
        t1 = -b/(2.0*a);
        t2 = INF;

    }else{
        
        t1 = (-b + sqrt(delta))/(2.0*a);
        t2 = (-b - sqrt(delta))/(2.0*a);
        if(t1 < -0.001){
            t1 = INF;
        }
        if(t2 < -0.001){
            t2 = INF;
        }
    }
    double taux = t1;
    t1 = (taux<t2)?taux:t2;
    t2 = (taux<t2)?t2:taux;
  
    inter.collision.x = e.x + t1 * d.x;
    inter.collision.y = e.y + t1 * d.y;
    inter.collision.z = e.z + t1 * d.z;
    POINT aux;
    aux.x = inter.collision.x - cylinder->center.x;
    aux.y = inter.collision.y - cylinder->center.y;
    aux.z = inter.collision.z - cylinder->center.z;

    double dist = aux.x * cylinder->axis.x + aux.y * cylinder->axis.y + aux.z * cylinder->axis.z;
    inter.m = dist / sqrt(myPow(cylinder->axis.x,2) + myPow(cylinder->axis.y,2) + myPow(cylinder->axis.z,2));
    if(inter.m >= cylinder->d1 && inter.m <= cylinder->d2){
        inter.t = t1;
        inter.object = obj;
    }else{
        inter.collision.x = e.x + t2 * d.x;
        inter.collision.y = e.y + t2 * d.y;
        inter.collision.z = e.z + t2 * d.z;
        aux.x = inter.collision.x - cylinder->center.x;
        aux.y = inter.collision.y - cylinder->center.y;
        aux.z = inter.collision.z - cylinder->center.z;

        dist = aux.x * cylinder->axis.x + aux.y * cylinder->axis.y + aux.z * cylinder->axis.z;
        inter.m = dist / sqrt(myPow(cylinder->axis.x,2) + myPow(cylinder->axis.y,2) + myPow(cylinder->axis.z,2));
        if(inter.m >= cylinder->d1 && inter.m <= cylinder->d2){
            inter.t = t2;
            inter.object = obj;
        }
    }
            
    
    return inter;
}

INTERSECTION IntersectionCone(void* obj, POINT e, POINT d){

    if(DEBUG)printf("IntersectionCylinder\n");
    CONE* cone = (CONE*)((OBJ*)obj)->object;
    
    INTERSECTION inter;

    POINT X;
    X.x = e.x - cone->center.x;
    X.y = e.y - cone->center.y;
    X.z = e.z - cone->center.z;

    inter.t == INF;
    double t1,t2;
    double a = myPow(d.x,2) + myPow(d.y,2) + myPow(d.z,2) - (1 + myPow(cone->angle,2)) * myPow((d.x*cone->axis.x + d.y*cone->axis.y + d.z*cone->axis.z),2);

    double b = 2.0 * (
        ((d.x )*(X.x) + (d.y )*(X.y) + (d.z)*(X.z)) - 
        (1 + myPow(cone->angle,2)) * 
        (d.x*cone->axis.x + d.y*cone->axis.y + d.z*cone->axis.z) * 
        ((cone->axis.x )*(X.x) + (cone->axis.y)*(X.y) + (cone->axis.z)*(X.z)) );

    double g = 
            myPow((X.x), 2.0) + myPow((X.y), 2.0) + myPow((X.z), 2.0) -
            (1 + myPow(cone->angle,2)) * 
            myPow(((cone->axis.x )*(X.x) + (cone->axis.y)*(X.y) + (cone->axis.z)*(X.z)), 2.0);

    double delta = myPow(b, 2.0) - 4.0 * g * a;

    //printf("%lf\n",delta );
    if(delta < -0.001){
        //printf("no hay Colision con Cilindro\n");
        t1 = INF;
        t2 = INF;
    }else if(delta < 0.001){
        t1 = -b/(2.0*a);
        t2 = INF;

    }else{
        
        t1 = (-b + sqrt(delta))/(2.0*a);
        t2 = (-b - sqrt(delta))/(2.0*a);
        if(t1 < -0.001){
            t1 = INF;
        }
        if(t2 < -0.001){
            t2 = INF;
        }
    }
    double taux = t1;
    t1 = (taux<t2)?taux:t2;
    t2 = (taux<t2)?t2:taux;
  
    inter.collision.x = e.x + t1 * d.x;
    inter.collision.y = e.y + t1 * d.y;
    inter.collision.z = e.z + t1 * d.z;
    POINT aux;
    aux.x = inter.collision.x - cone->center.x;
    aux.y = inter.collision.y - cone->center.y;
    aux.z = inter.collision.z - cone->center.z;

    double dist = aux.x * cone->axis.x + aux.y * cone->axis.y + aux.z * cone->axis.z;
    inter.m = dist / sqrt(myPow(cone->axis.x,2) + myPow(cone->axis.y,2) + myPow(cone->axis.z,2));
    if(inter.m >= cone->d1 && inter.m <= cone->d2){
        inter.t = t1;
        inter.object = obj;
    }else{
        inter.collision.x = e.x + t2 * d.x;
        inter.collision.y = e.y + t2 * d.y;
        inter.collision.z = e.z + t2 * d.z;
        aux.x = inter.collision.x - cone->center.x;
        aux.y = inter.collision.y - cone->center.y;
        aux.z = inter.collision.z - cone->center.z;

        dist = aux.x * cone->axis.x + aux.y * cone->axis.y + aux.z * cone->axis.z;
        inter.m = dist / sqrt(myPow(cone->axis.x,2) + myPow(cone->axis.y,2) + myPow(cone->axis.z,2));
        if(inter.m >= cone->d1 && inter.m <= cone->d2){
            inter.t = t2;
            inter.object = obj;
        }
    }
    
    
    return inter;
}

INTERSECTION IntersectionDisc(void* obj, POINT e, POINT d){
    if(DEBUG)printf("IntersectionSphere\n");
    DISC* disc = (DISC*)((OBJ*)obj)->object;
    INTERSECTION inter;
    long double D = -(disc->plane->normal.x * disc->plane->center.x + disc->plane->normal.y * disc->plane->center.y + disc->plane->normal.z * disc->plane->center.z);

    inter.t =
            ((disc->plane->normal.x * e.x + disc->plane->normal.y * e.y + disc->plane->normal.z * e.z) + D )
            / (disc->plane->normal.x * d.x + disc->plane->normal.y * d.y + disc->plane->normal.z * d.z);

    inter.t = -inter.t;
    if(inter.t < -0.001)inter.t = INF;
    else{

        inter.collision.x = e.x + inter.t * d.x;
        inter.collision.y = e.y + inter.t * d.y;
        inter.collision.z = e.z + inter.t * d.z;
        double dist = sqrt(myPow(inter.collision.x - disc->plane->center.x,2) + myPow(inter.collision.y - disc->plane->center.y,2) + myPow(inter.collision.z - disc->plane->center.z,2));
        //printf("Disco %lf\n", dist);
        if(dist < disc->radius1 ||  dist > disc->radius2)inter.t = INF;
        if(inter.t != INF){
            //printf("Disco %lf\n", inter.t);
            inter.object = obj;
        }
    }
    
    
    return inter;
}

INTERSECTION IntersectionOval(void* obj, POINT e, POINT d){
    if(DEBUG)printf("IntersectionSphere\n");
    OVAL* oval = (OVAL*)((OBJ*)obj)->object;
    INTERSECTION inter;
    long double D = -(oval->normal.x * oval->center1.x + oval->normal.y * oval->center1.y + oval->normal.z * oval->center1.z);

    inter.t =
            ((oval->normal.x * e.x + oval->normal.y * e.y + oval->normal.z * e.z) + D )
            / (oval->normal.x * d.x + oval->normal.y * d.y + oval->normal.z * d.z);

    inter.t = -inter.t;
    if(inter.t < -0.001)inter.t = INF;
    else{

        inter.collision.x = e.x + inter.t * d.x;
        inter.collision.y = e.y + inter.t * d.y;
        inter.collision.z = e.z + inter.t * d.z;
        double dist1 = sqrt(myPow(inter.collision.x - oval->center1.x,2) + myPow(inter.collision.y - oval->center1.y,2) + myPow(inter.collision.z - oval->center1.z,2));
        double dist2 = sqrt(myPow(inter.collision.x - oval->center2.x,2) + myPow(inter.collision.y - oval->center2.y,2) + myPow(inter.collision.z - oval->center2.z,2));
        double dist = dist1 + dist2;
        //printf("ovalo %lf\n", dist);
        if(dist > oval->radius)inter.t = INF;
        if(inter.t != INF){
            //printf("Disco %lf\n", inter.t);
            inter.object = obj;
        }
    }
    
    
    return inter;
}

POINT GetNormalSphere(INTERSECTION inter){
    
	SPHERE* sphere = (SPHERE*)((OBJ*)inter.object)->object;
	POINT N;
	N.x = inter.collision.x - sphere->center.x;
    N.y = inter.collision.y - sphere->center.y;
    N.z = inter.collision.z - sphere->center.z;
    double n = sqrt(myPow(N.x, 2) + myPow(N.y, 2) + myPow(N.z, 2));
    N.x /=n;
    N.y /=n;
    N.z /=n;
    return N;
}

POINT GetNormalPoly(INTERSECTION inter){
    
	POLYGON * polygon = (POLYGON *)((OBJ*)inter.object)->object;
	return polygon->plane->normal;
}

POINT GetNormalPlane(INTERSECTION inter){
    
	PLANE * plane = (PLANE *)((OBJ*)inter.object)->object;
	return plane->normal;
}

POINT GetNormalDisc(INTERSECTION inter){
    
    DISC * disc = (DISC *)((OBJ*)inter.object)->object;
    return disc->plane->normal;
}

POINT GetNormalOval(INTERSECTION inter){
    
    OVAL * oval = (OVAL *)((OBJ*)inter.object)->object;
    return oval->normal;
}

POINT GetNormalCylinder(INTERSECTION inter){
    
    CYLINDER* cylinder = (CYLINDER*)((OBJ*)inter.object)->object;
    //printf("%lf\n", inter.m);
    
    POINT N;
    N.x = inter.collision.x - cylinder->axis.x * inter.m - cylinder->center.x;
    N.y = inter.collision.y - cylinder->axis.y * inter.m - cylinder->center.y;
    N.z = inter.collision.z - cylinder->axis.z * inter.m - cylinder->center.z;
    double n = sqrt(myPow(N.x, 2) + myPow(N.y, 2) + myPow(N.z, 2));
    N.x /=n;
    N.y /=n;
    N.z /=n;
    return N;
}

POINT GetNormalCone(INTERSECTION inter){
    
    CONE* cone = (CONE*)((OBJ*)inter.object)->object;
    //printf("%lf\n", inter.m);
    
    POINT N;
    N.x = inter.collision.x - cone->center.x - (1 + myPow(cone->angle,2)) * cone->axis.x * inter.m;
    N.y = inter.collision.y - cone->center.y - (1 + myPow(cone->angle,2)) * cone->axis.y * inter.m;
    N.x = inter.collision.z - cone->center.z - (1 + myPow(cone->angle,2)) * cone->axis.z * inter.m;
    double n = sqrt(myPow(N.x, 2) + myPow(N.y, 2) + myPow(N.z, 2));
    N.x /=n;
    N.y /=n;
    N.z /=n;
    return N;
}

void loadScene(){
    scanf("H_SIZE %d\n",&H_SIZE);
    scanf("V_SIZE %d\n",&V_SIZE);
    scanf("SHOWPROGRESS %d\n",&SHOWPROGRESS);
    scanf("ANTIALIASING %d\n",&ANTIALIASING);
    scanf("SHADOWS %d\n",&SHADOWS);
    scanf("eye x %lf y %lf z %lf\n",&eye.x,&eye.y,&eye.z);
    scanf("N_LIGHTS %d\n",&N_LIGHTS);

    lights = malloc(sizeof(LIGHT*)*N_LIGHTS);

    for (int i = 0; i < N_LIGHTS; ++i)
    {
        lights[i] = malloc(sizeof(LIGHT));
        scanf("x %lf y %lf z %lf\n",&lights[i]->pos.x,&lights[i]->pos.y,&lights[i]->pos.z);
        scanf("intensity %lf\n",&lights[i]->intensity);
        lights[i]->pos.y *=-1;
    }
    POINT points[3];
    POINT A,B;
    scanf("N_OBJECTS %d\n",&N_OBJECTS);
    objects = malloc(sizeof(OBJ)*N_OBJECTS);
    int n_points;
    int type;
    for (int i = 0; i < N_OBJECTS; ++i)
    {
        scanf("type %d\n",&type);
        switch(type){
            case 1: // Spheres
                objects[i].object = (void*)malloc(sizeof(SPHERE));
                SPHERE* auxSphere = (SPHERE*)objects[i].object;
                scanf("radius %lf\n",&auxSphere->radius);
                scanf("x %lf y %lf z %lf\n",&auxSphere->center.x,&auxSphere->center.y,&auxSphere->center.z);
                scanf("R %lf G %lf B %lf\n",&objects[i].color.R,&objects[i].color.G,&objects[i].color.B);
                //printf("R %lf G %lf B %lf\n",objects[i].color.R,objects[i].color.G,objects[i].color.B);
                scanf("Ka %lf Kd %lf Ks %lf Kn %lf\n",&objects[i].Ka,&objects[i].Kd,&objects[i].Ks,&objects[i].Kn);
                objects[i].fun_ptr = &IntersectionSphere;
                objects[i].norm_ptr = &GetNormalSphere;
                break;
            case 2: // Planes
                objects[i].object = (void*)malloc(sizeof(PLANE));
                PLANE* auxPlane = (PLANE*)objects[i].object;
                

                for (int k = 0; k < 3; ++k)
                {
                scanf("x %lf y %lf z %lf\n", &points[k].x, &points[k].y, &points[k].z);
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
                scanf("R %lf G %lf B %lf\n",&objects[i].color.R,&objects[i].color.G,&objects[i].color.B);
                scanf("Ka %lf Kd %lf Ks %lf Kn %lf\n",&objects[i].Ka,&objects[i].Kd,&objects[i].Ks,&objects[i].Kn);
                break;
            case 3: // Polygons
                
                scanf("N_POINTS %d\n",&n_points);
               

                objects[i].object = (void*)malloc(sizeof(POLYGON));
                POLYGON* auxPolygon = (POLYGON*)objects[i].object;
                auxPolygon->points = malloc(sizeof(POINT)*n_points);
                auxPolygon->n_points = n_points;

                for (int k = 0; k < n_points; ++k)
                {
                    auxPolygon->points[k] = malloc(sizeof(POINT));
                    scanf("x %lf y %lf z %lf\n", &auxPolygon->points[k]->x, &auxPolygon->points[k]->y, &auxPolygon->points[k]->z);
                    

                }

                auxPolygon->plane = malloc(sizeof(PLANE));

                //Se crean los vectores de direccion
                if(n_points >= 3){

                    A.x = auxPolygon->points[1]->x - auxPolygon->points[0]->x;
                    A.y = auxPolygon->points[1]->y - auxPolygon->points[0]->y;
                    A.z = auxPolygon->points[1]->z - auxPolygon->points[0]->z;
             
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

                scanf("R %lf G %lf B %lf\n",&objects[i].color.R,&objects[i].color.G,&objects[i].color.B);
                scanf("Ka %lf Kd %lf Ks %lf Kn %lf\n",&objects[i].Ka,&objects[i].Kd,&objects[i].Ks,&objects[i].Kn);
                break;
            case 4: // Cylinders
                objects[i].object = (void*)malloc(sizeof(CYLINDER));
                CYLINDER* auxCylinder = (CYLINDER*)objects[i].object;

                objects[i].fun_ptr = &IntersectionCylinder;
                objects[i].norm_ptr = &GetNormalCylinder;

                scanf("radius %lf\n",&auxCylinder->radius);
                scanf("d1 %lf\n",&auxCylinder->d1);
                scanf("d2 %lf\n",&auxCylinder->d2);
                scanf("center x %lf y %lf z %lf\n",&auxCylinder->center.x,&auxCylinder->center.y,&auxCylinder->center.z);
                scanf("axis x %lf y %lf z %lf\n",&auxCylinder->axis.x,&auxCylinder->axis.y,&auxCylinder->axis.z);

                scanf("R %lf G %lf B %lf\n",&objects[i].color.R,&objects[i].color.G,&objects[i].color.B);
                scanf("Ka %lf Kd %lf Ks %lf Kn %lf\n",&objects[i].Ka,&objects[i].Kd,&objects[i].Ks,&objects[i].Kn);
                break;
            case 5: // Disc
                objects[i].object = (void*)malloc(sizeof(DISC));
                DISC* auxDisc = (DISC*)objects[i].object;
                auxDisc->plane = malloc(sizeof(PLANE));
                
                scanf("radius1 %lf\n",&auxDisc->radius1);
                scanf("radius2 %lf\n",&auxDisc->radius2);
                scanf("center x %lf y %lf z %lf\n",&auxDisc->plane->center.x,&auxDisc->plane->center.y,&auxDisc->plane->center.z);
                for (int k = 0; k < 3; ++k)
                {
                scanf("x %lf y %lf z %lf\n", &points[k].x, &points[k].y, &points[k].z);
                }

                //Se crean los vectores de direccion
               
                A.x = points[1].x - points[0].x;
                A.y = points[1].y - points[0].y;
                A.z = points[1].z - points[0].z;
              
                B.x = points[2].x - points[0].x;
                B.y = points[2].y - points[0].y;
                B.z = points[2].z - points[0].z;

                objects[i].fun_ptr = &IntersectionDisc;
                objects[i].norm_ptr = &GetNormalDisc;
                auxDisc->plane->normal = dot(A, B);
                auxDisc->plane->normal.x *= -1.0;
                auxDisc->plane->normal.y *= -1.0;
                auxDisc->plane->normal.z *= -1.0;
                scanf("R %lf G %lf B %lf\n",&objects[i].color.R,&objects[i].color.G,&objects[i].color.B);
                scanf("Ka %lf Kd %lf Ks %lf Kn %lf\n",&objects[i].Ka,&objects[i].Kd,&objects[i].Ks,&objects[i].Kn);
                break;
            case 6: // Oval
                objects[i].object = (void*)malloc(sizeof(OVAL));
                OVAL* auxOval = (OVAL*)objects[i].object;
                scanf("radius %lf\n",&auxOval->radius);
                scanf("center1 x %lf y %lf z %lf\n",&auxOval->center1.x,&auxOval->center1.y,&auxOval->center1.z);
                scanf("center2 x %lf y %lf z %lf\n",&auxOval->center2.x,&auxOval->center2.y,&auxOval->center2.z);
                for (int k = 0; k < 3; ++k)
                {
                scanf("x %lf y %lf z %lf\n", &points[k].x, &points[k].y, &points[k].z);
                }

                //Se crean los vectores de direccion
               
                A.x = points[1].x - points[0].x;
                A.y = points[1].y - points[0].y;
                A.z = points[1].z - points[0].z;
              
                B.x = points[2].x - points[0].x;
                B.y = points[2].y - points[0].y;
                B.z = points[2].z - points[0].z;

                objects[i].fun_ptr = &IntersectionOval;
                objects[i].norm_ptr = &GetNormalOval;
                auxOval->normal = dot(A, B);
                auxOval->normal.x *= -1.0;
                auxOval->normal.y *= -1.0;
                auxOval->normal.z *= -1.0;
                scanf("R %lf G %lf B %lf\n",&objects[i].color.R,&objects[i].color.G,&objects[i].color.B);
                scanf("Ka %lf Kd %lf Ks %lf Kn %lf\n",&objects[i].Ka,&objects[i].Kd,&objects[i].Ks,&objects[i].Kn);
                break;
            case 7: // Cone
                objects[i].object = (void*)malloc(sizeof(CONE));
                CONE* auxCone = (CONE*)objects[i].object;

                objects[i].fun_ptr = &IntersectionCone;
                objects[i].norm_ptr = &GetNormalCone;

                scanf("angle %lf\n",&auxCone->angle);
                scanf("d1 %lf\n",&auxCone->d1);
                scanf("d2 %lf\n",&auxCone->d2);
                scanf("center x %lf y %lf z %lf\n",&auxCone->center.x,&auxCone->center.y,&auxCone->center.z);
                scanf("axis x %lf y %lf z %lf\n",&auxCone->axis.x,&auxCone->axis.y,&auxCone->axis.z);

                scanf("R %lf G %lf B %lf\n",&objects[i].color.R,&objects[i].color.G,&objects[i].color.B);
                scanf("Ka %lf Kd %lf Ks %lf Kn %lf\n",&objects[i].Ka,&objects[i].Kd,&objects[i].Ks,&objects[i].Kn);
                break;
            default:
                break;
        }
    }
}

void init() {
    srand (time(NULL));

    AmbientIlluminationIntensity = 0.2;

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
    double tmin;
    INTERSECTION intersection;
    intersection.t = INF;
    INTERSECTION auxIntersection;

    for (int i = 0; i < N_OBJECTS; ++i)
    {
        if(DEBUG)printf("------First_Intersection %d\n", i);

        auxIntersection = (objects[i].fun_ptr)((void*)&objects[i], e, d);

        if(auxIntersection.t > (double)SHADOW_K && auxIntersection.t < intersection.t){
            intersection = auxIntersection;
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
		//printf("interseccion R %lf G %lf B %lf\n",color.R,color.G,color.B);
        //POINT intersectionPoint;
        double intensity = 0;
        double E = 0;

        intersection.collision.x = e.x + intersection.t * d.x;
        intersection.collision.y = e.y + intersection.t * d.y;
        intersection.collision.z = e.z + intersection.t * d.z;
        double n;
        POINT N = (obj->norm_ptr)(intersection);
        POINT L;

        for (int i = 0; i < N_LIGHTS; ++i)
        {


    		double Fatt,cosNL,cosVR;

            L.x = lights[i]->pos.x - intersection.collision.x;
            L.y = lights[i]->pos.y - intersection.collision.y;
            L.z = lights[i]->pos.z - intersection.collision.z;
            n = sqrt(myPow(L.x, 2) + myPow(L.y, 2) + myPow(L.z, 2));
            L.x /=n;
            L.y /=n;
            L.z /=n;


            bool lid = true;
            if(SHADOWS){
            	INTERSECTION intersectionLight;
            	intersectionLight = First_Intersection(intersection.collision, L);

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
            		if(cosNL > 0)intensity += (cosNL * obj->Kd * lights[i]->intensity * Fatt);
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
    	if(SHOWPROGRESS)printf("%lf\n", ((double)i / (double)H_SIZE)*100);
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

	            double L = sqrt(myPow(w.x - eye.x, 2) + myPow(w.y - eye.y, 2) + myPow(w.z - eye.z, 2));
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

			            double L = sqrt(myPow(w.x - eye.x, 2) + myPow(w.y - eye.y, 2) + myPow(w.z - eye.z, 2));
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
