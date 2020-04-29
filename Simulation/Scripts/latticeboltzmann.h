#ifndef LATTICEBOLTZMANN_H
#define LATTICEBOLTZMANN_H

#include <iostream>
#include <fstream>
#include <cmath>

/*
  Para proportion=12 se manejaron valores de lambda=20 y A=150.
  Además de que se comienza a graficar en x=5, para no ver la fuente inicial.
*/


const int Lr=20;
const int Lt=4;
const int Lz=50;

const int Q=19;
const double W0=1.0/3;
//Constante de reflexión
const double rmax= 200;
const double rmin=2;
const double zmin=0;
const double zmax= 1200;
const double thetamax=2*M_PI;
const double thetamin=0;
const double lambda= 4*zmax/9.0;

const double C=0.5; // C<0.707 celdas/click
const double C2=C*C;
const double TresC2=3*C*C;
const double AUX0=1-TresC2*(1-W0);
const double Cs=0.5;
const double Cs2=Cs*Cs;
const double UCs2=1.0/Cs2;
const double D=0.65;
const double tau=1;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
const double deltar= (rmax-rmin)/Lr;
const double Udeltar=1.0/deltar;
const double deltatheta= (thetamax-thetamin)/Lt;
const double Udeltat=1.0/deltatheta;
const double deltaz=(zmax-zmin)/Lz;
const double Udeltaz=1.0/deltaz;
const double deltas=deltar*deltatheta*deltaz;
class LatticeBoltzmann
{
private:
  double w[Q];
  int V[3][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy, V[2][i]=V_iz
  double f[Lr][Lt][Lz][Q], fnew[Lr][Lt][Lz][Q]; // f[ix][iy][iz][iz][i]
public:
  LatticeBoltzmann(void);
  double r(int i);	
  double theta(int j);
  double z(int k);
  double Mg(int i, int j);
  double Sig(int i , int j, int k);
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
  double GetRho(int ix, int iy, int iz, bool algo);
};

#endif // LATTICEBOLTZMANN_H
