#include "latticeboltzmann.h"
#include "omp.h"

LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0] = W0; w[1] = w[2] = w[3] = w[4] = (1 - W0) / 4.0;;
  //Cargar los vectores
  V[0][0] = 0; V[0][1] = 1;  V[0][2] = 0;  V[0][3] = -1;  V[0][4]=  0;
  V[1][0] = 0; V[1][1] = 0;  V[1][2] = 1;  V[1][3] = 0;  V[1][4]= -1;
}

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  int i; double suma = 0;
    for(i=0;i<Q;i++){
    if(UseNew) suma += fnew[ix][iy][i]; else suma += f[ix][iy][i];
    }
  return suma;
}

double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew){
  int i; double suma;
    for(suma=0,i=0;i<Q;i++){
      if(UseNew) suma += fnew[ix][iy][i]*V[0][i]; else suma += f[ix][iy][i]*V[0][i];
    }
    return suma;
}

double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma += fnew[ix][iy][i]*V[1][i]; else suma += f[ix][iy][i]*V[1][i];
  return suma;
}


double LatticeBoltzmann::feq(double rho0, double Jx0, double Jy0, int i){
  if(i==0)
    return rho0 * AUX0;
  else
    return w[i] * (TresC2 * rho0 + 3* (V[0][i] *Jx0 + V[1][i] * Jy0));
}

void LatticeBoltzmann::Colisione(void){
  int ix,iy,iz,i; double rho0,Jx0,Jy0;  //Para cada celda

#pragma omp paralel for
  {
  for(iy=0;iy<Ly;iy++){
    for(ix=0;ix<Lx;ix++){
        //Calcular las cantidades macroscópicas
	rho0 = rho(ix, iy, false);  Jx0 = Jx(ix, iy, false);  Jy0 = Jy(ix, iy, false);  
        for(i=0; i<Q; i++){
	 fnew[ix][iy][i] = UmUtau*f[ix][iy][i] + Utau*feq(rho0, Jx0, Jy0, i);
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
      for(int i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i] = fnew[ix][iy][i];
}
}
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  #pragma omp paralel for
  {
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
        f[ix][iy][i]=feq(rho0, Jx0, Jy0, i);
}
}
void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix,iy; double lambda,omega,rho0,Jx0,Jy0;
  lambda = 10; omega = 2*M_PI / lambda; ix = Lx/2; iy = Ly/2;
  rho0 = 10 * sin(omega*t); Jx0 = Jx(ix, iy, false); Jy0 = Jy(ix, iy, false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][i] = feq(rho0, Jx0, Jy0, i);
}

void LatticeBoltzmann::Imprimase(const char * NombreArchivo, double t){
  std::ofstream MiArchivo(NombreArchivo);
  double rho0,Jx0,Jy0;
  #pragma omp paralel for
  {
  MiArchivo<< t<<'\n';
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
        rho0 = rho(ix, iy, false);  //Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
        MiArchivo<< ix << " " << iy << " " <<rho0<<'\n';
      //MiArchivo<<'\n';
    }
    //MiArchivo<<'\n';
  }
  MiArchivo<<'\n';
  }
  if (t==999){MiArchivo.close();}
}

void LatticeBoltzmann::Imprimir(int t, int ix, int iy, const char * NombreArchivo){
  double rho0 = rho(ix, iy, false);
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << '\t' << rho0 << '\n';
  ofs.close();
}



