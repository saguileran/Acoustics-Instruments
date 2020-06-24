#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "omp.h"


const int proportion = 1;
const int Lx = 501*proportion, Ly = 50*proportion;
const int LFx = 330*(proportion), LFy = 14*(proportion);

const double ke = 0, k2 = 0, kF = 1; //k_enviorment = 0 means wall absortion
const double Aperture_x = 2*proportion;
const double Hole_pos = LFx/3;

const int Q = 5;
const double W0 = 1.0 / 3;
const double k = 1; //Constante de reflexión

const double C = 0.5; // C<0.707 celdas/click
const double TresC2 = 3 * C * C;
const double AUX0 = 1 - TresC2 * (1 - W0);

const double tau = 0.5;
const double Utau = 1.0 / tau;
const double UmUtau = 1 - Utau;


//----------- CLASS ------------------------
class LatticeBoltzmann
{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, bool UseNew);
  double  Jx(int ix, int iy, bool UseNew);
  double  Jy(int ix, int iy, bool UseNew);
  double  Jz(int ix, int iy, bool UseNew);
  double feq(double rho0, double Jx0, double Jy0, int i);
  void Colisione(void);
  void Adveccione(void);
  void Inicie(double rho0, double Jx0, double Jy0);
  void ImponerCampos(int t);
  void Imprimase(const char * NombreArchivo, int t);
  void Imprimir(int t, int ix, int iy, const char * NombreArchivo);
  void Microphone(int t, int ix, int iy, const char * NombreArchivo);
};


//---------------------FUNCTIONS------------------------------
LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0] = W0;   w[1] = w[2] = w[3] = w[4] = (1 - W0) / 4.0;
  //Cargar los vectores
  V[0][0] = 0; V[0][1] = 1;  V[0][2] = 0;  V[0][3] = -1;  V[0][4]=  0;
  V[1][0] = 0; V[1][1] = 0;  V[1][2] = 1;  V[1][3] =  0;  V[1][4]= -1;
}

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  int i; double suma = 0;
    for(i=0;i<Q;i++){if(UseNew) suma += fnew[ix][iy][i]; else suma += f[ix][iy][i]; }
  return suma;
}

double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew){
  int i; double suma;
    for(suma=0,i=0;i<Q;i++){if(UseNew) suma += fnew[ix][iy][i]*V[0][i]; else suma += f[ix][iy][i]*V[0][i]; }
    return suma;
}

double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++){if(UseNew) suma += fnew[ix][iy][i]*V[1][i]; else suma += f[ix][iy][i]*V[1][i];}
  return suma;
}


double LatticeBoltzmann::feq(double rho0, double Jx0, double Jy0, int i){
  if(i==0){ return rho0 * AUX0;}
  else{     return w[i] * (TresC2 * rho0 + 3* (V[0][i] *Jx0 + V[1][i] * Jy0));}
}

void LatticeBoltzmann::Colisione(void){
  int ix, iy, iz, i; double rho0, Jx0, Jy0;  //for all cell

#pragma omp paralel for
  {
  for(iy=0;iy<Ly;iy++){
    for(ix=0;ix<Lx;ix++){
        //Calcular las cantidades macroscópicas
	rho0 = rho(ix, iy, false);  Jx0 = Jx(ix, iy, false);  Jy0 = Jy(ix, iy, false);
	
	fnew[ix][iy][0] = UmUtau*f[ix][iy][0] + Utau*feq(rho0, Jx0, Jy0, 0);
	//for(i=0; i<Q; i++){ fnew[ix][iy][i] = UmUtau*f[ix][iy][i] + Utau*feq(rho0, Jx0, Jy0, i);}

	//Left flute wall
	if(ix == 20 &&  iy >= Ly/2 - LFy/2 && iy <= Ly/2 + LFy/2 ) {fnew[ix][iy][1] = kF * f[ix][iy][2]; fnew[ix][iy][2] = kF * f[ix][iy][1];}
	else if(ix == Lx - 1 || ix == 1){ fnew[ix][iy][1] = ke * fnew[ix][iy][3]; fnew[ix][iy][3] = ke * fnew[ix][iy][1]; } 
	else{ fnew[ix][iy][1] = UmUtau*f[ix][iy][1] + Utau*feq(rho0, Jx0, Jy0, 1);
    	      fnew[ix][iy][3] = UmUtau*f[ix][iy][3] + Utau*feq(rho0, Jx0, Jy0, 3); }

	//hirizontall fulte walls, the comment lines are the holes
	if(((iy == Ly/2 - LFy/2 && ix >= 20 && ix <= 20 + LFx) || (iy == Ly/2 + LFy/2  &&  ix >= 20 && ix <= 20 + LFx)
	    //&&  not(ix >= 20 + Hole_pos + LFx/3 - Aperture_x/2 && ix <= 20 + Hole_pos + LFx/3 + Aperture_x/2 && iy == Ly/2 + LFy/2)
	    //&&  not(ix >= 20 + Hole_pos - Aperture_x/2 && ix <= 20 + Hole_pos + Aperture_x/2 && iy == Ly/2 + LFy/2)
	    )
	   )
	  {fnew[ix][iy][4] =  kF * fnew[ix][iy][2]; fnew[ix][iy][4] =  kF * fnew[ix][iy][4];}
	else if(iy == Ly - 2 || iy == 1){ fnew[ix][iy][2] = ke *  f[ix][iy][4]; fnew[ix][iy][4] = ke *  f[ix][iy][2]; } 
	else{ fnew[ix][iy][2] = UmUtau*f[ix][iy][2] + Utau*feq(rho0, Jx0, Jy0, 2);
	      fnew[ix][iy][4] = UmUtau*f[ix][iy][4] + Utau*feq(rho0, Jx0, Jy0, 4); }
	
	//std::cout << ix << " " << iy << " " << rho0 << std::endl; //microphone place
	//if(ix == Lx-1 and iy == Ly-1){std::cout << " " << std::endl;}
	
    }
  }
 }
}
  
void LatticeBoltzmann::Adveccione(void){
  #pragma omp paralel for
  {
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)	
	//This condition disable lattice periodicity
	if( ix + V[0][i] < Lx && ix + V[0][i] >= 0 && iy + V[1][i] < Ly && iy + V[1][i] >= 0){ 
	  f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i] = fnew[ix][iy][i];}
  }
}
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  #pragma omp paralel for
  {
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[ix][iy][i] = feq(rho0, Jx0, Jy0, i);
}
}
void LatticeBoltzmann::ImponerCampos(int t){
  int i, ix, iy; double lambda, omega, rho0, Jx0, Jy0;
  
  //sin(omega * t), declare initial function variables
  lambda = 10; omega = 2 * M_PI / lambda; ix = 22 ; iy = (Ly) / 2;

  //Initialice macroscopic cuantities
  rho0 = 10; //* sin(omega*t);
  Jx0 = Jx(ix, iy, false);  Jy0 = Jy(ix, iy, false);

  //make a pulse of 500 steps
  if(t < 2000){
    //for(iy=(Ly-LFy+2)/2;iy<(Ly+LFy)/2;iy++){
    for(i=0;i<Q;i++)
      fnew[ix][iy][i] = feq(rho0, Jx0, Jy0, i);
    //    }
  }
}

void LatticeBoltzmann::Imprimase(const char * NombreArchivo, int t){
  std::ofstream MiArchivo(NombreArchivo + std::to_string(t));
  double rho0, Jx0, Jy0;
  #pragma omp paralel for
  {
  MiArchivo << "X,Y,rho\n";
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
        rho0 = rho(ix, iy, false);  //Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
	//export data to paraview visualization
	MiArchivo << ix << "," << iy << "," << rho0 << '\n';
    }
  }
  }
  MiArchivo.close();
}

void LatticeBoltzmann::Imprimir(int t, int ix, int iy, const char * NombreArchivo){
  double rho0 = rho(ix, iy, false);
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << '\t' << rho0 << '\n';
  ofs.close();
}

void LatticeBoltzmann::Microphone(int t, int ix, int iy, const char * NombreArchivo){
  double suma = 0; 

  for(int i = 1; i < 4; i++){
    double rho0 = rho(ix + i, iy - 1, false);
    double rho1 = rho(ix + i, iy,     false);
    double rho2 = rho(ix + i, iy + 1, false);

    suma += rho0 + rho1 + rho2;
    }
  
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << '\t' << suma / 9 << '\n';
  ofs.close();
}




//--------MAIN---------------------------------
int main(void){
  
  LatticeBoltzmann Ondas;  
  int t,tmax=10000;
  
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
  */
  
  Ondas.Inicie(0,0,0);
  for(t=0;t<tmax;t++){

    //print actual step
    std::cout << t << std::endl;

    //Making lattice steps
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();

    //Export microphoes data, time vs pressure
    Ondas.Imprimir(t, 20+LFx,    Ly/2, "LongPulse10k-0mm.dat");
    Ondas.Imprimir(t, 20+LFx+10, Ly/2, "LongPulse10k-10mm.dat");
    Ondas.Imprimir(t, 20+LFx+60, Ly/2, "LongPulse10k-60mm.dat");


    //Commands to make data animation
    if(t%5 == 0){Ondas.Imprimase("LongPulse10k.csv.", t);}

    //Uncomment to use GNUPLOT animation
    //std::cout << "splot 'Ondas.dat' using 1:2:3  with points palette pointsize 3 pointtype 7 " << std::endl;
  }
  return 0;
}
