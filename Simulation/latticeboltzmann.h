#ifndef LATTICEBOLTZMANN_H
#define LATTICEBOLTZMANN_H

#include <iostream>
#include <fstream>
#include <cmath>

/*
  Para proportion=12 se manejaron valores de lambda=20 y A=150.
  Además de que se comienza a graficar en x=5, para no ver la fuente inicial.
*/


const int proportion=3;
const int Lx=20;
const int Ly=4;
const int Lz=100;

const int Q=19;
const double W0=1.0/3;
//Constante de reflexión
const double rmax= 200;
const double rmin=2;
const double zmin=0;
const double zmax= 1200;
const double lambda= 4*zmax/9.0;

const double C=0.5; // C<0.707 celdas/click
const double C2=C*C;
const double TresC2=3*C*C;
const double AUX0=1-TresC2*(1-W0);
const double Cs=0.5;
const double Cs2=Cs*Cs;
const double D=0.65;
const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
const double deltar= (rmax-rmin)/Lx;
const double deltatheta= 2*M_PI/Ly;
const double deltatheta2=deltatheta*deltatheta;
const double deltaz=zmax/Lz;
const double deltas=deltaz*deltar*deltatheta;
class LatticeBoltzmann
{
private:
  double w[Q];
  int V[3][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy, V[2][i]=V_iz
  double f[Lx][Ly][Lz][Q], fnew[Lx][Ly][Lz][Q]; // f[ix][iy][iz][iz][i]
public:
  LatticeBoltzmann(void);
  double r(int i);	
  double rho(int ix,int iy,int iz,bool UseNew);
  double Jx(int ix,int iy,int iz,bool UseNew);
  double Jy(int ix,int iy,int iz,bool UseNew);
  double Jz(int ix,int iy,int iz,bool UseNew);
  double feq(double rho0,double Jx0,double Jy0,double Jz0,int i,int ix);
  void Colisione(void);
  void Adveccione(void);
  void Inicie(double rho0,double Jx0,double Jy0, double Jz0);
  void ImponerCampos(int t);
  void Imprimase(const char * NombreArchivo);
  void Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo);
  bool Columna(int x1, int x2, int x);
  bool Rectangulo(int x1, int x2, int x, int y1, int y2, int y);
  double GetRho(int ix, int iy, int iz, bool algo){double rho0 = rho(ix,iy,iz,algo);  return rho0;};
};

#endif // LATTICEBOLTZMANN_H
