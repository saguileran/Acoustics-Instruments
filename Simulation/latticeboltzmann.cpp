#include <omp.h>
#include "latticeboltzmann.h"

LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=w[5]=w[6]=1.0/18; w[7]=w[8]=w[9]=w[10]=w[11]=w[12]=w[13]=w[14]=w[15]=w[16]=w[17]=w[18]=1.0/36;
  //Cargar los vectores
 
  V[0][0]=0;  V[1][0]=0;  V[2][0]=0;

  V[0][1]=1;  V[0][2]=-1;  V[0][3]=0;  V[0][4]=0;  V[0][5]=0;  V[0][6]=0; 
  V[1][1]=0;  V[1][2]=0;  V[1][3]=1;  V[1][4]=-1; V[1][5]=0;  V[1][6]=0;
  V[2][1]=0;  V[2][2]=0;  V[2][3]=0;  V[2][4]=0;  V[2][5]=1;  V[2][6]=-1;

  V[0][7]=1; V[0][8]=-1; V[0][9]=1; V[0][10]=-1; V[0][11]=1; 
  V[1][7]=1; V[1][8]=-1; V[1][9]=-1; V[1][10]=1; V[1][11]=0;  
  V[2][7]=0; V[2][8]=0;  V[2][9]=0; V[2][10]=0;  V[2][11]=1; 

  V[0][12]=-1; V[0][13]=1;  V[0][14]=-1; V[0][15]=0; V[0][16]=0; V[0][17]=0;
  V[1][12]=0;  V[1][13]=0;  V[1][14]=0; V[1][15]=1; V[0][16]=-1; V[1][17]=1;
  V[2][12]=-1; V[2][13]=-1; V[2][14]=1; V[2][15]=1; V[0][16]=-1; V[2][17]=-1;

  V[0][18]=0; V[1][18]=-1; V[2][18]=1;
}
double LatticeBoltzmann::r(int i)
{
  double r=i*deltar+rmin;
  return r;
}
double LatticeBoltzmann::rho(int ix,int iy, int iz,bool UseNew){
  int i; double suma,r0; r0=r(ix);
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][iz][i]; else suma+=f[ix][iy][iz][i];
  suma=suma/(r0*deltas);
  return suma;
}

double LatticeBoltzmann::Jx(int ix,int iy, int iz,bool UseNew){
  int i; double suma,suma2,r0, r0,suma3;r0=r(ix);
  for(suma=0,i=0;i<Q;i++)
    if(UseNew){ suma+=fnew[ix][iy][iz][i]*V[0][i];  suma2+=fnew[ix][iy][iz][i];suma3=w[i]*V[0][i];} else{ suma+=f[ix][iy][iz][i]*V[0][i];suma2+=f[ix][iy][iz][i];suma3=w[i]*V[0][i];}
  suma+=0.5*suma2*C2/(deltar*r0);
  suma+=0.5*(1/Cs2)*suma3*deltas*r0*(Cs2-(C2/(deltar*deltar)));
  suma=suma/(r0*deltas);
  return suma;
}

double LatticeBoltzmann::Jy(int ix,int iy,int iz,bool UseNew){
   int i; double suma,r0, r0,suma3;r0=r(ix);
  for(suma=0,i=0;i<Q;i++)
    if(UseNew){ suma+=fnew[ix][iy][iz][i]*V[1][i];suma3=w[i]*V[1][i];} else{ suma+=f[ix][iy][iz][i]*V[1][i];suma3=w[i]*V[1][i];}
  suma+=0.5*(1/Cs2)*suma3*deltas*r0*(Cs2-(C2/(deltatheta*deltatheta*r0*r0)));
  suma=suma/(r0*deltas);
  return suma;
}

double LatticeBoltzmann::Jz(int ix,int iy,int iz,bool UseNew){
  int i; double suma,r0, r0,suma3;r0=r(ix);
  for(suma=0,i=0;i<Q;i++)
    if(UseNew){ suma+=fnew[ix][iy][iz][i]*V[2][i];suma3=w[i]*V[2][i];} else{ suma+=f[ix][iy][iz][i]*V[2][i];suma3=w[i]*V[2][i];}
  suma+=0.5*(1/Cs2)*suma3*deltas*r0*(Cs2-(C2/(deltaz*deltaz)));
  suma=suma/(r0*deltas);
  return suma;
}

double LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,double Jz0,int i,int ix){
  double r0=r(ix);
  if(i==0)
    return w[i]*rho0*deltas*r0;
  else
    return deltas*r0*w[i]*(rho0+(V[0][i]*Jx0+V[1][i]*Jy0+V[2][i]Jz0)/Cs2);
}

void LatticeBoltzmann::Colisione(void){
  int ix,iy,i; double rho0,Jx0,Jy0;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
      //Calcular las cantidades macroscópicas

	rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false);
	if((ix*ix)+(iy*iy)<(rmax*rmax)){ 
	  for (i=0;i<Q;i++)
	    fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,i,ix);}
	else{
       

          fnew[ix][iy][iz][2]=D*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=D*f[ix][iy][iz][2];
          fnew[ix][iy][iz][4]=D*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=D*f[ix][iy][iz][4];
          fnew[ix][iy][iz][6]=D*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=D*f[ix][iy][iz][6];
	}
    }
}
void LatticeBoltzmann::Adveccione(void){
 for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][(iz+V[2][i]+Lz)%Ly][i]=fnew[ix][iy][iz][i];
}

void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0, double Jz0){
  int ix,iy,iz,i;
#pragma omp parallel for private(ix, iy, iz, i)
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
        for(i=0;i<Q;i++)
          f[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
}

void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix,iy,iz,A; double lambda,omega,rho0,Jx0,Jy0,Jz0;
  omega=2*M_PI*C/lambda; A=2;
  ix=0; iy=0; iz=0;
  
  rho0=A*sin(omega*t)*rmin; Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
}

void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  std::ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0,Jz0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      //for(int iz=0;iz<Lz;iz++){
        int iz=20;
        rho0=rho(ix,iy,iz,true);  //Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
        MiArchivo<<ix<<" "<<iy<<" "<<rho0<<'\n';
      //}
      //MiArchivo<<'\n';
    }
    MiArchivo<<'\n';
  }
  MiArchivo.close();
}

void LatticeBoltzmann::Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo){
  double rho0 = rho(ix,iy,iz,true);
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << '\t' << rho0 << '\n';
  ofs.close();
}

bool LatticeBoltzmann::Columna(int x1, int x2, int x)
{
  if(x>=x1 && x<=x2){return true;}
  else{return false;}
}

bool LatticeBoltzmann::Rectangulo(int x1, int x2, int x, int y1, int y2, int y)
{
  if(x>=x1 && x<=x2){
    if(y>=y1 && y<=y2){return true;}
    else{return false;}
  }
  else{return false;}
}
