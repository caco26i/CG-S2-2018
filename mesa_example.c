/*
 * Instituto Tecnologico de Costa Rica
 * Escuela de Ingenieria en Computacion
 * Computer Graphics
 *
 * Programa: Mesa Example
 * Archivo:  mesa_example.c
 */

#include "mesa_example.h"
#include "malloc.h"

void draw_scene ();

COLOR **buffer;

int main(int argc, char** argv) 
{
  int i, j;

  buffer = (COLOR **)malloc(H_SIZE * sizeof(COLOR*));
  for (i = 0; i < H_SIZE; i++) 
      {
       buffer[i] = (COLOR *)malloc(V_SIZE * sizeof(COLOR));
      }

  for (i = 0; i < H_SIZE; i++) 
      {
       for (j = 0; j < V_SIZE; j++) 
           {
            buffer[i][j].r = 0;
            buffer[i][j].g = 0;
            buffer[i][j].b = 0;
           }
      }

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(H_SIZE,V_SIZE);
  glutCreateWindow("Mesa Example");
  glClear(GL_COLOR_BUFFER_BIT);
  gluOrtho2D(-0.5, H_SIZE +0.5, -0.5, V_SIZE + 0.5);
  glutDisplayFunc(draw_scene);
  glutMainLoop();
}


void draw_scene() {
  static int last_x = 0;
  int i, j;
  COLOR color;

  for (i = 0; i < last_x; i++) 
      {
       for (j = 0; j < V_SIZE; j++) 
           {
            glColor3f (buffer[i][j].r,buffer[i][j].g,buffer[i][j].b);
            glBegin (GL_POINTS);
            glVertex2i (i,j);
            glEnd();
           }
      }

  for (i = last_x; i < H_SIZE; i++) 
      {
       for (j = 0; j < V_SIZE; j++) 
         {
          buffer[i][j].r = (double)(i % (H_SIZE / 10)) / (double)(H_SIZE / 10);
          buffer[i][j].g = (double)(j % (V_SIZE / 10)) / (double)(V_SIZE / 10);
          buffer[i][j].b = (double)(i) / (double)(H_SIZE);
          glColor3f (buffer[i][j].r,buffer[i][j].g,buffer[i][j].b);
          glBegin(GL_POINTS);
          glVertex2i(i,j);
          glEnd();
          last_x = i;
         }
      }

  glFlush();
}




