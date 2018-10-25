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

void loadScene(){
	scanf("H_SIZE %d\n",&H_SIZE);
	scanf("V_SIZE %d\n",&V_SIZE);
	scanf("ANTIALIASING %d\n",&ANTIALIASING);
	scanf("SHADOWS %d\n",&SHADOWS);
	scanf("N_SPHERES %d\n",&N_SPHERES);

	printf("H_SIZE %d\n",H_SIZE);
	printf("V_SIZE %d\n",V_SIZE);
	printf("ANTIALIASING %d\n",ANTIALIASING);
	printf("SHADOWS %d\n",SHADOWS);
	printf("N_SPHERES %d\n",N_SPHERES);

	spheres = malloc(sizeof(SPHERE*)*N_SPHERES);

	for (int i = 0; i < N_SPHERES; ++i)
	{
		spheres[i] = malloc(sizeof(SPHERE));
	    scanf("radius %f\n",&spheres[i]->radius);
	    scanf("x %f\n",&spheres[i]->center.x);
	    scanf("y %f\n",&spheres[i]->center.y);
	    scanf("z %f\n",&spheres[i]->center.z);
	    scanf("R %f\n",&spheres[i]->color.R);
	    scanf("G %f\n",&spheres[i]->color.G);
	    scanf("B %f\n",&spheres[i]->color.B);
	    scanf("Ka %f\n",&spheres[i]->Ka);
	    scanf("Kd %f\n",&spheres[i]->Kd);
	    scanf("Ks %f\n",&spheres[i]->Ks);
	    scanf("Kn %f\n",&spheres[i]->Kn);

	    printf("radius %f\n",spheres[i]->radius);
	    printf("x %f\n",spheres[i]->center.x);
	    printf("y %f\n",spheres[i]->center.y);
	    printf("z %f\n",spheres[i]->center.z);
	    printf("R %f\n",spheres[i]->color.R);
	    printf("G %f\n",spheres[i]->color.G);
	    printf("B %f\n",spheres[i]->color.B);
	    printf("Ka %f\n",spheres[i]->Ka);
	    printf("Kd %f\n",spheres[i]->Kd);
	    printf("Ks %f\n",spheres[i]->Ks);
	    printf("Kn %f\n",spheres[i]->Kn);
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

    background.R = 0.03;
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
    float b = 2.0 * ((d.x )*( e.x - sphere->center.x) + (d.y )*( e.y - sphere->center.y) +
                     (d.z )*( e.z - sphere->center.z));
    float g = myPow((e.x - sphere->center.x), 2.0) + myPow((e.y - sphere->center.y), 2.0) +
              myPow((e.z - sphere->center.z), 2.0) - myPow(sphere->radius, 2.0);
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
    intersection.object = (void*)sphere;
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
        	if(auxIntersection.t >= 0 && auxIntersection.t < intersection.t)intersection = auxIntersection;
        }
        
    }
    return intersection;
}


COLOR De_que_color(POINT e, POINT d) {
    COLOR color;
    INTERSECTION intersection;
    intersection = First_Intersection(e, d);

    if (intersection.t == INF){
        color = background;
    }

    else {
        SPHERE *obj = (SPHERE *) intersection.object;
        color = obj->color;

        POINT intersectionPoint;
        float intensity = 0;
        float E = 0;

        intersectionPoint.x = e.x + intersection.t * d.x;
        intersectionPoint.y = e.y + intersection.t * d.y;
        intersectionPoint.z = e.z + intersection.t * d.z;
        float n;
        POINT N;
        N.x = intersectionPoint.x - obj->center.x;
        N.y = intersectionPoint.y - obj->center.y;
        N.z = intersectionPoint.z - obj->center.z;
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
            	
	            if (intersectionLight.t != INF){
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

				

				if(cosVR > 0)E += (myPow(cosVR,obj->Kn) * obj->Ks * lights[i]->intensity * Fatt);
	            if(cosNL > 0)intensity += (cosNL * obj->Kd * lights[i]->intensity * Fatt);
            }
            
            //printf("Fatt %f\n", n );
            
            
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

    printf("Fin del programa %s...\n\n", argv[0]);


    return 1;
}
