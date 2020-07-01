#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "omp.h"


const int proportion = 1;
const int Lx = 40*proportion, Ly = 40*proportion, Lz = 40 * proportion;

//const int LFx = 330*(proportion), LFy = 14*(proportion);

const double k1 = 0, k2 = 0, kF = 0;

const int Q = 7, dt = 1;
const double W0 = 1.0 / 4;
const double k = 1; //Constante de reflexión

const double C = 0.5, Cs = 1; // C<0.707 celdas/click
const double C1 = (C / Cs) * (C / Cs);
const double C2 = 1 + C1 * (W0 - 1);

const double tau  = 0.5;
const double tau1 = dt / tau;
const double tau2 = 1 - tau1;


//----------- CLASS ------------------------
class LatticeBoltzmann
{
private:
  double w[Q];
  int V[3][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Lz][Q], fnew[Lx][Ly][Lz][Q]; // f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, int iz, bool UseNew);
  double Jx(int ix, int iy, int iz, bool UseNew);
  double Jy(int ix, int iy, int iz, bool UseNew);
  double Jz(int ix, int iy, int iz, bool UseNew);
  double feq(double rho0, double Jx0, double Jy0, double Jz0, int i);
  void Colisione(void);
  void Adveccione(void);
  void Inicie(double rho0, double Jx0, double Jy0, double Jz0);
  void ImponerCampos(int t);
  void Imprimase(const char * NombreArchivo, int t);
  void Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo);
};


//---------------------FUNCTIONS------------------------------
LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0] = W0;    w[1] = w[2] = w[3] = w[4] = w[5] = w[6] = W0 / 2.0;
  //Cargar los vectores
  V[0][0] = 0; V[0][1] = 1;  V[0][2] = -1;  V[0][3] = 0;   V[0][4]=  0;  V[0][5]=  0;  V[0][6]=  0;
  V[1][0] = 0; V[1][1] = 0;  V[1][2] =  0;  V[1][3] = 1;   V[1][4]= -1;  V[1][5]=  0;  V[1][6]=  0;
  V[2][0] = 0; V[2][1] = 0;  V[2][2] =  0;  V[2][3] = 0;   V[2][4]=  0;  V[2][5]=  1;  V[2][6]=  -1;
}

double LatticeBoltzmann::rho(int ix, int iy, int iz, bool UseNew){
  int i; double suma = 0;
  for(i=0;i<Q;i++){
    if(UseNew) suma += fnew[ix][iy][iz][i]; else suma += f[ix][iy][iz][i];
  }
  return suma;
}

double LatticeBoltzmann::Jx(int ix, int iy, int iz, bool UseNew){
  int i; double suma;
    for(suma=0,i=0;i<Q;i++){
      if(UseNew) suma += fnew[ix][iy][iz][i]*V[0][i]; else suma += f[ix][iy][iz][i]*V[0][i];
    }
    return suma;
}

double LatticeBoltzmann::Jy(int ix, int iy, int iz, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++){
    if(UseNew) suma += fnew[ix][iy][iz][i]*V[1][i]; else suma += f[ix][iy][iz][i]*V[1][i];
  }
  return suma;
}

double LatticeBoltzmann::Jz(int ix, int iy, int iz, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++){ if(UseNew) suma += fnew[ix][iy][iz][i]*V[2][i]; else suma += f[ix][iy][iz][i]*V[2][i];}
  return suma;
}


double LatticeBoltzmann::feq(double rho0, double Jx0, double Jy0, double Jz0, int i){
  if(i==0)
    return rho0 * C2;
  else
    return w[i] * C1 * rho0 + (V[0][i] *Jx0 + V[1][i] * Jy0 + V[2][i] * Jz0);
}

void LatticeBoltzmann::Colisione(void){
  int ix,iy,iz,i; double rho0,Jx0,Jy0, Jz0;  //Para cada celda

#pragma omp paralel for
  {
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iz<Lz;iz++){
      for(iz=0;iz<Lz;iz++){
        //Calcular las cantidades macroscópicas
	rho0 = rho(ix, iy, iz, false);  Jx0 = Jx(ix, iy, iz, false);  Jy0 = Jy(ix, iy, iz, false); Jz0 = Jy(ix, iy, iz, false);
	//fnew[ix][iy][0] = tau2*f[ix][iy][iz][0] + tau1*feq(rho0, Jx0, Jy0, Jz0, 0);
	for(i=0; i<Q; i++){
	  fnew[ix][iy][iz][i] = tau2 * f[ix][iy][iz][i] + tau1 * feq(rho0, Jx0, Jy0, Jz0, i);
	}
      }
    }
  }
  }
}
  
void LatticeBoltzmann::Adveccione(void){
  #pragma omp paralel for
  {
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	for(int i=0;i<Q;i++)
	  if( ix + V[0][i] < Lx && ix + V[0][i] >= 0 && iy + V[1][i] < Ly && iy + V[1][i] >= 0 && iz + V[2][i] < Lz && iz + V[2][i] >= 0){ //This condition disable lattice periodicity
	    f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][(iz+V[2][i]+Lz)%Lz][i] = fnew[ix][iy][iz][i];}
  }
}
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0, double Jz0){
  #pragma omp paralel for
  {
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	for(int i=0;i<Q;i++)
	  f[ix][iy][iz][i] = feq(rho0, Jx0, Jy0, Jz0, i);
  }
}
void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix,iy, iz; double lambda,omega,rho0,Jx0,Jy0, Jz0;
  lambda = 10; omega = 2*M_PI / lambda; ix = Lx/2 ; iy = Ly/2; iz = Lz/2;
  rho0 = 10 * sin(omega*t);
  Jx0 = Jx(ix, iy, iz, false); Jy0 = Jy(ix, iy, iz, false); Jz0 = Jz(ix, iy, iz, false);
  //  if(t < 1000){
  //  for(iy=(Ly-LFy+2)/2;iy<(Ly+LFy)/2;iy++){
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i] = feq(rho0, Jx0, Jy0, Jz0, i);
      //}
  //}
}

void LatticeBoltzmann::Imprimase(const char * NombreArchivo, int t){
  std::ofstream MiArchivo(NombreArchivo  + std::to_string(t) + ".csv");
  double rho0, Jx0, Jy0, Jz0;
  MiArchivo <<  "X,Y,Z,rho" << '\n';
  #pragma omp paralel for
  {
  MiArchivo << '\n';
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      for(int iz=0;iz<Lz;iz++){
        rho0 = rho(ix, iy, iz, false);  //Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
        MiArchivo << ix << "," << iy << "," << iz << "," << rho0 << '\n';
      }
    }
  }
  }
  MiArchivo.close();
}

void LatticeBoltzmann::Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo){
  double rho0 = rho(ix, iy, iz, false);
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << '\t' << rho0 << '\n';
  ofs.close();
}




//--------MAIN---------------------------------
int main(void){
  
  LatticeBoltzmann Ondas;  
  int t,tmax=100;
  
  //GNUPLOT
  /*
  // Estos comandos se descomentan si se quiere guardar el gif
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set output 'pelicula0.gif'" << std::endl;
  
  //Estos comandos se descomentan para hacer el gif
  std::cout << "set pm3d map; set palette color positive" << std::endl;
  std::cout << "set palette defined (-1 \"red\", 0 \"white\", 1 \"blue\")" << std::endl;
  std::cout << "set cbrange[-1:1]" << std::endl;
  std::cout << "set xrange[-1:501]; set yrange[-1:51]; set zrange[-1:1]" << std::endl;
  //std::cout << "set view map scale 1 " << std::endl;

  //std::cout << "set view map;  set size ratio .9 " << std::endl;
  //std::cout << "set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back " << std::endl;
  //std::cout << "set object 1 rect fc rgb 'black' fillstyle solid 1.0 " << std::endl;
  //std::cout << " " << std:endl;
  //std::cout << " " << std:endl;
  */
  
  Ondas.Inicie(0,0,0,0);
  for(t=0;t<tmax;t++){
    std::cout << t << std::endl;
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
    //Ondas.Imprimir(t,25,25,25,"datos.dat");     //Este comando se tiene para graficar la amplitud en función del tiempo en el punto x,y,z

    //Estos comandos son los que permiten hacer el gif
    Ondas.Imprimase("Ondas", t); 
    //std::cout << "splot 'Ondas.dat' using 1:2:3  with points palette pointsize 3 pointtype 7 " << std::endl;
  }
  //Ondas.Imprimase("Ondas.dat");

  return 0;
}
