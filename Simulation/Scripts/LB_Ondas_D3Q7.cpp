#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/gl.h>
#include <GL/glut.h>
#include "latticeboltzmann.h"

double rho1[Lr][Lt][Lz];
LatticeBoltzmann Ondas;
int t,t_otro;



int main(int argc, char** argv)
{
 
  //Gnuplot

  int t,tmax=200;
  Ondas.Inicie(0,0,0,0);
  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    if(t<=(10)){
      Ondas.ImponerCampos(t);
    }
    Ondas.Adveccione();
    std::cout << t << '\t' ;
    std::cout << Ondas.GetRho(rmax/10,0,0,true) << '\t'; //Primera posición
    std::cout << '\n';
  }
  
  return 0;
}
