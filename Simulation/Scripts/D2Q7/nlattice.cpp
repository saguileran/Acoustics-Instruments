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
  glPointSize(3.0);
  glBegin(GL_POINTS);
   #pragma omp paralel for
  {
  for(int ix=0; ix<Lx; ix++){
    for(int iy=0; iy<Ly; iy++){
	rho0 = rho1[ix][iy];
	glColor3f(0.0, 1.0-rho0*20.0, rho0*20.0);
	glVertex3f(ix*0.01,iy*0.01,0.01);
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
  glEnd();
  std::cout << t << std::endl;
  glFlush();
}

void init(void)
{
  // select clearing (background) color
  glClearColor(0.0, 0.0, 0.0, 0.0);

  //initialize viewing values */
  GLfloat aspect = (GLfloat)width / (GLfloat)height;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glTranslatef(-0.25f, -0.1f, -0.0f);  // Move into the screen
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
  
  //Gnuplot

  int t,tmax=1000;
  /*
  // Estos comandos se descomentan si se quiere guardar el gif
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set output 'pelicula0.gif'" << std::endl;
  //Estos comandos se descomentan para hacer el gif
  std::cout << "set pm3d" << std::endl;
  std::cout << "set palette defined (-1 \"red\", 0 \"white\", 1 \"blue\")" << std::endl;
  std::cout << "set cbrange[-1:1]" << std::endl;
  std::cout << "set xrange[-1:41]; set yrange[-1:41]; set zrange[-1:5]" << std::endl;
  */
  //Ondas.Inicie(0,0,0,0);
  /*
  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
    //Este comando se tiene para graficar la amplitud en función del tiempo en el punto x,y,z
    //Ondas.Imprimir(t,25,25,25,"Datos.dat");
    //Estos comandos son los que permiten hacer el gif
    
    Ondas.Imprimase("Ondas1.dat", t);
    //std::cout << "splot 'Ondas.dat'" << std::endl;
    Ondas.Imprimir*(
  }
  //  Ondas.Imprimase("Datos.dat");
  */
  return 0;
}
