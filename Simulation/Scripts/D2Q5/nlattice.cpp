#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/gl.h>
#include <GL/glut.h>
#include "latticeboltzmann.h"
#include "omp.h"

double rho1[Lx][Ly];
LatticeBoltzmann Ondas;
int t,t_otro;
int height= 1400, width = 1400;

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT);

  double rho0;  //int iz=Lz-proportion;
  Ondas.Colisione();
  Ondas.ImponerCampos(t);
  Ondas.Imprimase("Ondas1.dat", t);
  Ondas.Adveccione();
  // Ondas.Imprimir(t, Lx/2, Ly/2, "OndasData.dat");

  glEnd();
  std::cout << t << " " << Ondas.rho(20 + LFx, Ly/2, false) << std::endl;
  glFlush();
  
  glPointSize(3.0);
  glBegin(GL_POINTS);
   #pragma omp paralel for
  {
  for(int ix=0; ix<Lx; ix++){
    for(int iy=0; iy<Ly; iy++){
	rho0 = rho1[ix][iy];
	glColor3f(0.0, 1.0-rho0*20.0, rho0*20.0);
	if(((ix == 20 &&  iy >= Ly/2 - LFy/2 && iy <= Ly/2 + LFy/2) || (iy == Ly/2 - LFy/2 &&  ix >= 20 && ix <= 20 + LFx) || (iy == Ly/2 + LFy/2  &&  ix >= 20 && ix <= 20 + LFx))
	  // &&  not(ix >= 20 + Hole_pos - Aperture_x/2 && ix <= 20 + Hole_pos + Aperture_x/2 && iy == Ly/2 + LFy/2)
	   )
	  { glVertex3f(0, 0, 0);}
	else{
	  glVertex3f(ix*0.01,iy*0.01,0.01);
	}
	/* Setup the view of the cube. */
	glMatrixMode(GL_PROJECTION);
	gluPerspective( /* field of view in degree */ 120.0,
			/* aspect ratio */ 2.0,
	  			/* Z near */ 1.0, /* Z far */ 10.0);
	//  glMatrixMode(GL_MODELVIEW);
	//      gluLookAt(20.0, 50.0, 50.0,  /* eye is at (0,0,5) */
	// 10.0, 10.0, 10.0,      /* center is at (0,0,0) */
	//		  0.0, 1.0, 0.);      /* up is in positive Y direction */
	
    }
  }
}
  //glEnd();
  //std::cout << t << std::endl;
  //glFlush();
}

void init(void)
{
  // select clearing (background) color
  glClearColor(0.0, 0.0, 0.0, 0.0);

  //initialize viewing values */
  GLfloat aspect = (GLfloat)width / (GLfloat)height;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glTranslatef(-1.0f, -1.0f, -0.0f);  // Move into the screen
  gluLookAt(0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 1.0, 0.0);
}

void AmplitudDisplay(void)
{
  #pragma omp paralel for
  {
  for(int ix=0; ix<Lx; ix++){
    for(int iy=0; iy<Ly; iy++){
        rho1[ix][iy] = Ondas.rho(ix, iy, true);
    }
  }
  t++;
  glutPostRedisplay();
}
}
void mouse(int button, int state, int x, int y)
{
  switch (button) {
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN)
      glutIdleFunc(AmplitudDisplay);
      break;
    case GLUT_MIDDLE_BUTTON:
      if (state == GLUT_DOWN)
      glutIdleFunc(NULL);
      break;
    default:
      break;
  }
}


int main(int argc, char** argv)
{
  //OpenGL

  Ondas.Inicie(0,0,0);
  
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(width, height);
  glutInitWindowPosition(700, 700);
  glutCreateWindow("Animation");
  init();
  glutDisplayFunc(display);
  glutMouseFunc(mouse);
  glutMainLoop();
  
  return 0;
}
