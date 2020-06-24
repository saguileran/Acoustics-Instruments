#ifndef LATTICEBOLTZMANN_H
#define LATTICEBOLTZMANN_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "omp.h"
/*
  Para proportion=12 se manejaron valores de lambda=20 y A=150.
  Además de que se comienza a graficar en x=5, para no ver la fuente inicial.
*/
const int proportion = 1;
const int Lx = 34*proportion, Ly = 2*proportion, Lz = 2*proportion;
const int LFx = 32*(proportion), LFy = 16*(proportion);

const double k1 = 0, k2 = 0, kF = 0;
const double Aperture_x = 2*proportion;
const double Hole_pos = LFx/3;

const int Q = 7;
//const double W0 = 1.0 / 3;
const double k = 1; //Constante de reflexión

const double C = 0.5, Cs = 1; // C<0.707 celdas/click
const double dt = 1 ; //simulation time
const double C1 = 1 + (C**2 / cs**2) * (W[0] - 1);
const double C2 = w[i] / Cs**2;

const double tau = 0.5;
const double tau1 = 1 - dt / tau;
const double tau2 = dt / tau;

class LatticeBoltzmann
{
private:
  double w[Q];
  int V[3][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy, V[2][i]=V_iy
  double f[Lx][Ly][Lz][Q], fnew[Lx][Ly][Lz][Q]; // f[ix][iy][iz][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, int iz, bool UseNew);
  double Jx(int ix, int iy, int iz, bool UseNew);
  double Jy(int ix, int iy, int iz, bool UseNew);
  double Jz(int ix, int iy, int iz,  bool UseNew);
  double feq(double rho0, double Jx0, double Jy0, double Jz0, int i);
  void Colisione(void);
  void Adveccione(void);
  void Inicie(double rho0, double Jx0, double Jy0, double Jz0);
  void ImponerCampos(int t);
  void Imprimase(const char * NombreArchivo, double t);
  void Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo);
};

#endif // LATTICEBOLTZMANN_H
