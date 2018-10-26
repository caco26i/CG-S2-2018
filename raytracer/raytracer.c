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

int H_SIZE;
int V_SIZE;
int N_POLYGONS;
int N_SPHERES;
int N_LIGHTS;
int ANTIALIASING;
int SHADOWS;

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
    dot.x = a.z*b.y - a.y*b.z;
    dot.y = -(a.x*b.z - a.z*b.z);
    dot.z = a.x*b.y - a.y*b.x;
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

    int cases;

    if (fabs (poly->plane->normal.x) > fabs (poly->plane->normal.y) && fabs (poly->plane->normal.x) > fabs (poly->plane->normal.z) ){
        cases = 0;
    } else if (fabs (poly->plane->normal.y) > fabs (poly->plane->normal.x) && fabs (poly->plane->normal.y) > fabs (poly->plane->normal.z) ){
        cases = 1;
    } else {
        cases = 2;
    }

//    testx = punto.x;
//    testy = punto.y;

    for (int i = 0; i < poly->n_points; ++i){
        if (cases == 0){
            points[i].x = poly->points[i]->y - punto.y;
            points[i].y = poly->points[i]->z - punto.z;
            testx = punto.y;
            testy = punto.z;
        } else if (cases == 1) {
            points[i].x = poly->points[i]->x - punto.x;
            points[i].y = poly->points[i]->z - punto.z;
            testx = punto.x;
            testy = punto.z;
        }
        else {
            points[i].x = poly->points[i]->x - punto.x;
            points[i].y = poly->points[i]->y - punto.y;
            testx = punto.x;
            testy = punto.y;
        }
    }

    int i, j, c = 0;
    for (i = 0, j = poly->n_points-1; i < poly->n_points; j = i++) {
        if ( ((points[i].y>testy) != (points[j].y>testy)) &&
             (testx < (points[j].x-points[i].x) * (testy-points[i].y) / (points[j].y-points[i].y) + points[i].x) )
            c = !c;
    }
    return c;
}

void loadScene(){
	scanf("H_SIZE %d\n",&H_SIZE);
	scanf("V_SIZE %d\n",&V_SIZE);
	scanf("ANTIALIASING %d\n",&ANTIALIASING);
	scanf("SHADOWS %d\n",&SHADOWS);

	scanf("N_POLYGONS %d\n",&N_POLYGONS);
	printf("H_SIZE %d\n",H_SIZE);
	printf("V_SIZE %d\n",V_SIZE);
	printf("ANTIALIASING %d\n",ANTIALIASING);
	printf("SHADOWS %d\n",SHADOWS);
	printf("N_POLYGONS %d\n",N_POLYGONS);

	polygons = malloc(sizeof(POLYGON*)*N_POLYGONS);


    for (int i = 0; i < N_POLYGONS; ++i) {
        POINT a;
        POINT b;
        POINT c;
        int n_points;
        scanf("N_POINTS %d\n",&n_points);
        printf("N_POINTS %d\n",n_points);

        polygons[i] = malloc(sizeof(POLYGON));
        polygons[i]->points = malloc(sizeof(POINT)*n_points);

        for (int k = 0; k < n_points; ++k)
        {
            polygons[i]->points[k] = malloc(sizeof(POINT));
            scanf("x %f y %f z %f\n", &polygons[i]->points[k]->x, &polygons[i]->points[k]->y, &polygons[i]->points[k]->z);
            printf("x %f y %f z %f\n", polygons[i]->points[k]->x, polygons[i]->points[k]->y, polygons[i]->points[k]->z);
        }

        polygons[i]->plane = malloc(sizeof(PLANE));
        //Se crean los vectores de direccion
        POINT A;
        A.x = polygons[i]->points[1]->x - polygons[i]->points[0]->x;
        A.y = polygons[i]->points[1]->y - polygons[i]->points[0]->y;
        A.z = polygons[i]->points[1]->z - polygons[i]->points[0]->z;
        POINT B;
        B.x = polygons[i]->points[2]->x - polygons[i]->points[0]->x;
        B.y = polygons[i]->points[2]->y - polygons[i]->points[0]->y;
        B.z = polygons[i]->points[2]->z - polygons[i]->points[0]->z;

        polygons[i]->plane->normal = dot(A, B);
        polygons[i]->plane->object.center.x = polygons[i]->points[0]->x;
        polygons[i]->plane->object.center.y = polygons[i]->points[0]->y;
        polygons[i]->plane->object.center.z = polygons[i]->points[0]->z;

        scanf("R %f\n",&polygons[i]->plane->object.color.R);
        scanf("G %f\n",&polygons[i]->plane->object.color.G);
        scanf("B %f\n",&polygons[i]->plane->object.color.B);
        scanf("Ka %f\n",&polygons[i]->plane->object.Ka);
        scanf("Kd %f\n",&polygons[i]->plane->object.Kd);
        scanf("Ks %f\n",&polygons[i]->plane->object.Ks);
        scanf("Kn %f\n",&polygons[i]->plane->object.Kn);

        printf("x %f\n",polygons[i]->plane->object.center.x);
        printf("y %f\n",polygons[i]->plane->object.center.y);
        printf("z %f\n",polygons[i]->plane->object.center.z);
        printf("x %f\n",polygons[i]->plane->normal.x);
        printf("y %f\n",polygons[i]->plane->normal.y);
        printf("z %f\n",polygons[i]->plane->normal.z);
        printf("R %f\n",polygons[i]->plane->object.color.R);
        printf("G %f\n",polygons[i]->plane->object.color.G);
        printf("B %f\n",polygons[i]->plane->object.color.B);
        printf("Ka %f\n",polygons[i]->plane->object.Ka);
        printf("Kd %f\n",polygons[i]->plane->object.Kd);
        printf("Ks %f\n",polygons[i]->plane->object.Ks);
        printf("Kn %f\n",polygons[i]->plane->object.Kn);
    }

    scanf("N_SPHERES %d\n",&N_SPHERES);
    spheres = malloc(sizeof(SPHERE*)*N_SPHERES);
    printf("N_SPHERES %d\n",N_SPHERES);

    for (int i = 0; i < N_SPHERES; ++i)
	{
		spheres[i] = malloc(sizeof(SPHERE));
	    scanf("radius %f\n",&spheres[i]->radius);
	    scanf("x %f\n",&spheres[i]->object.center.x);
	    scanf("y %f\n",&spheres[i]->object.center.y);
	    scanf("z %f\n",&spheres[i]->object.center.z);
	    scanf("R %f\n",&spheres[i]->object.color.R);
	    scanf("G %f\n",&spheres[i]->object.color.G);
	    scanf("B %f\n",&spheres[i]->object.color.B);
	    scanf("Ka %f\n",&spheres[i]->object.Ka);
	    scanf("Kd %f\n",&spheres[i]->object.Kd);
	    scanf("Ks %f\n",&spheres[i]->object.Ks);
	    scanf("Kn %f\n",&spheres[i]->object.Kn);

	    printf("radius %f\n",spheres[i]->radius);
	    printf("x %f\n",spheres[i]->object.center.x);
	    printf("y %f\n",spheres[i]->object.center.y);
	    printf("z %f\n",spheres[i]->object.center.z);
	    printf("R %f\n",spheres[i]->object.color.R);
	    printf("G %f\n",spheres[i]->object.color.G);
	    printf("B %f\n",spheres[i]->object.color.B);
	    printf("Ka %f\n",spheres[i]->object.Ka);
	    printf("Kd %f\n",spheres[i]->object.Kd);
	    printf("Ks %f\n",spheres[i]->object.Ks);
	    printf("Kn %f\n",spheres[i]->object.Kn);
	}

	scanf("N_LIGHTS %d\n",&N_LIGHTS);

	printf("N_LIGHTS %d\n",N_LIGHTS);

	lights = malloc(sizeof(LIGHT*)*N_LIGHTS);

	for (int i = 0; i < N_LIGHTS; ++i)
	{
		lights[i] = malloc(sizeof(LIGHT));
	    scanf("x %f\n",&lights[i]->pos.x);
	    scanf("y %f\n",&lights[i]->pos.y);
	    scanf("z %f\n",&lights[i]->pos.z);
	    scanf("intensity %f\n",&lights[i]->intensity);

	    printf("x %f\n",lights[i]->pos.x);
	    printf("y %f\n",lights[i]->pos.y);
	    printf("z %f\n",lights[i]->pos.z);
	    printf("intensity %f\n",lights[i]->intensity);
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


INTERSECTION IntersectionSphere(SPHERE *sphere, POINT e, POINT d){
    INTERSECTION intersection;
    //float a = pow((d.x - e.x), 2.0) + pow((d.y - e.y), 2.0) + pow((d.z - e.z), 2.0);
    float a = 1;
    float b = 2.0 * ((d.x )*( e.x - sphere->object.center.x) + (d.y )*( e.y - sphere->object.center.y) +
                     (d.z )*( e.z - sphere->object.center.z));
    float g = myPow((e.x - sphere->object.center.x), 2.0) + myPow((e.y - sphere->object.center.y), 2.0) +
              myPow((e.z - sphere->object.center.z), 2.0) - myPow(sphere->radius, 2.0);
    float delta = myPow(b, 2.0) - 4.0 * g * a;


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
    intersection.object = sphere->object;
    return intersection;
}

INTERSECTION IntersectionPoly(POLYGON *polygon, POINT e, POINT d){
    INTERSECTION intersection;
    PLANE* plane = polygon->plane;

    long double D = -(plane->normal.x * plane->object.center.x + plane->normal.y * plane->object.center.y + plane->normal.z * plane->object.center.z);

    intersection.t =
            ((plane->normal.x * e.x + plane->normal.y * e.y + plane->normal.z * e.z) + D )
            / (plane->normal.x * d.x + plane->normal.y * d.y + plane->normal.z * d.z);

    POINT intersection_point;
    intersection_point.x = e.x + intersection.t*d.x;
    intersection_point.y = e.y + intersection.t*d.y;
    intersection_point.z = e.z + intersection.t*d.z;

    if(intersection.t < -0.001 || esta_punto_poly(polygon, intersection_point))intersection.t = INF;
    intersection.object = plane->object;
    return intersection;
}

INTERSECTION First_Intersection(POINT e, POINT d) {
    float tmin;
    INTERSECTION intersection;
    intersection.t = INF;
    INTERSECTION auxIntersection;

    for (int i = 0; i < N_SPHERES; ++i)
    {
        auxIntersection = IntersectionSphere(spheres[i], e, d);
        if(auxIntersection.t > (float)SHADOW_K){
        	if(auxIntersection.t >= 0 && auxIntersection.t < intersection.t) intersection = auxIntersection;
        }

    }

    for (int i = 0; i < N_POLYGONS; ++i)
    {
        auxIntersection = IntersectionPoly(polygons[i], e, d);
        if(auxIntersection.t > (float)SHADOW_K){
            if(auxIntersection.t >= 0 && auxIntersection.t < intersection.t) intersection = auxIntersection;
        }

    }

    return intersection;
}

COLOR De_que_color(POINT e, POINT d) {
    INTERSECTION intersection;
    intersection = First_Intersection(e, d);
    COLOR color = intersection.object.color;

    if (intersection.t == INF){
        color = background;
    }

    else {
        OBJ obj = intersection.object;

        POINT intersectionPoint;
        float intensity = 0;
        float E = 0;

        intersectionPoint.x = e.x + intersection.t * d.x;
        intersectionPoint.y = e.y + intersection.t * d.y;
        intersectionPoint.z = e.z + intersection.t * d.z;
        float n;
        POINT N;
        N.x = intersectionPoint.x - obj.center.x;
        N.y = intersectionPoint.y - obj.center.y;
        N.z = intersectionPoint.z - obj.center.z;
        n = sqrt(myPow(N.x, 2) + myPow(N.y, 2) + myPow(N.z, 2));
        N.x /=n;
        N.y /=n;
        N.z /=n;

        POINT L;

        for (int i = 0; i < N_LIGHTS; ++i)
        {

        	INTERSECTION intersectionLight;
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
            	intersectionLight = First_Intersection(intersectionPoint, L);
            	
	            if (intersectionLight.t != INF && intersectionLight.t < n){
	            	lid = false;
			    }
            }

            if(lid){

            	
            	POINT R;
            	POINT V;

            	Fatt = fmin(1.0,1.0 / (myPow(n/100.0,2)));

            	cosNL = (L.x * N.x + L.y * N.y + L.z * N.z);

            	R.x = 2 * N.x * cosNL - L.x;
            	R.y = 2 * N.y * cosNL - L.y;
            	R.z = 2 * N.z * cosNL - L.z;
            	
            	V.x = -d.x;
            	V.y = -d.y;
            	V.z = -d.z;

            	cosVR = (V.x * R.x + V.y * R.y + V.z * R.z);

				
            	if(cosNL > 0){
            		if(cosVR > 0)E += (myPow(cosVR,obj.Kn) * obj.Ks * lights[i]->intensity * Fatt);
            		intensity += (cosNL * obj.Kd * lights[i]->intensity * Fatt);
            	}


            }
            
            //printf("Fatt %f\n", n );
            
            
        }

        intensity += obj.Ka * AmbientIlluminationIntensity;
        
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
    for (int i = 0; i < H_SIZE; i++) {
    	if(PROGRESS)printf("%f\n", ((float)i / (float)H_SIZE)*100);
        for (int j = 0; j < V_SIZE; j++) {
        	

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

