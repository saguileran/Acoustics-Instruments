#include <omp.h>
#include "latticeboltzmann.h"




//Pesos y velocidades




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

//Valores relacionados con la forma, coordenadas métrica etc



// coordenadas

double LatticeBoltzmann::r(int i)
{
  double r=i*deltar+rmin;
  return r;
}


double LatticeBoltzmann::theta(int j)
{
  double theta =j*deltatheta;
  return theta;
}

double LatticeBoltzmann::z(int k){
  double z= zmin-k*deltaz;
  return z;
}

//Tensor de métrica inverso g^{mu nu} 


double LatticeBoltzmann::Mg(int i, int j){
  double value;
  if (i==j){
    if (i==0)
      value=Udeltar*Udeltar;
    else if (i==1)
      value=r(i)*Udeltat*Udeltat;
    else
      value=Udeltaz*Udeltaz;
  }
  else value=0;
  return value;
}

// Tensor Christoffel symbols ke no se ke jeso


double LatticeBoltzmann::Sig(int i, int j, int k)
{
  if (i==j){
    if(i==1)
      if(k==0)
	return -deltatheta*deltatheta*Udeltar*r(k);}
  else if(k==j){
    if(j==1)
      if(i==0)
	return deltar/r(i);}   
  else if (k==i)
    { if (i==1)
	if(j==0)
	  return deltar/r(j);
    }
  else return 0;
}



// rho nuevo sin embargo no le quitamos el raiz de g sino que se hace el cambio en GetRho




double LatticeBoltzmann::rho(int ix,int iy, int iz,bool UseNew){
  int i; double suma,r0; r0=r(ix);
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][iz][i]; else suma+=f[ix][iy][iz][i];
  return suma;
}



double LatticeBoltzmann::GetRho(int ix, int iy, int iz, bool algo)
{
  double rho0 = rho(ix,iy,iz,algo); rho0=rho0/deltas; rho0=rho0/r(ix); return rho0;
}

// Densidades de flujo  J

double LatticeBoltzmann::Jx(int ix,int iy, int iz,bool UseNew){
  int i;
  double suma,raizg, r0,P; P=rho(ix,iy,iz,UseNew);
  r0=r(ix);
  raizg=deltas*r0;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew){
      suma+=fnew[ix][iy][iz][i]*V[0][i];
      suma+=0.5*w[i]*V[0][i]*raizg;
      for(int k =0; k<3; k++){
	suma-=w[i]*V[k][i]*Mg(0,k)*C2*UCs2;}
    }
  for( i=0; i<3; i++)
    for(int j=0; j<3;j++)
      {
	suma-=0.5*C2*Sig(i,j,0)*Mg(i,j);
      }
  return suma;
}

double LatticeBoltzmann::Jy(int ix,int iy,int iz,bool UseNew){
 int i;
  double suma,raizg, r0,P; P=rho(ix,iy,iz,UseNew);
  r0=r(ix);
  raizg=deltas*r0;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew){
      suma+=fnew[ix][iy][iz][i]*V[1][i];
      suma+=0.5*w[i]*V[1][i]*raizg;
      for(int k =0; k<3; k++){
	suma-=w[i]*V[k][i]*Mg(1,k)*C2*UCs2;}
    }
  for( i=0; i<3; i++)
    for(int j=0; j<3;j++)
      {
	suma-=0.5*C2*Sig(i,j,1)*Mg(i,j);
      }
  return suma;
}

double LatticeBoltzmann::Jz(int ix,int iy,int iz,bool UseNew){
  int i;
  double suma,raizg, r0,P; P=rho(ix,iy,iz,UseNew);
  r0=r(ix);
  raizg=deltas*r0;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew){
      suma+=fnew[ix][iy][iz][i]*V[0][i];
      suma+=0.5*w[i]*V[2][i]*raizg;
      for(int k =0; k<3; k++){
	suma-=w[i]*V[k][i]*Mg(2,k)*C2*UCs2;}
    }
  for( i=0; i<3; i++)
    for(int j=0; j<3;j++)
      {
	suma-=0.5*C2*Sig(i,j,2)*Mg(i,j);
      }
  return suma;
}

//funcion de equilibrio


double LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,double Jz0,int i,int ix){
  if(i==0)
    return w[i]*rho0;
  else
    return w[i]*(rho0+(V[0][i]*Jx0+V[1][i]*Jy0+V[2][i]*Jz0)*UCs2);
}

// Colisione, aqui hay mucha incertidumbre jaja si algo sale mal probablemente es de aqui

void LatticeBoltzmann::Colisione(void){
  int ix,iy,iz,i; double rho0,Jx0,Jy0,Jz0,r0,theta0,z0;
  //Para cada celda
  for(ix=0;ix<Lr;ix++)
    for(iy=0;iy<Lt;iy++)
      for(iz=0;iz<Lz;iz++){
      //Calcular las cantidades macroscópicas

	rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false);Jz0=Jz(ix,iy,iz,false);
	r0=r(ix); theta0=theta(iy); z0=z(iz);

	  for (i=0;i<Q;i++){
	    fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,Jz0,i,ix);}

      }}


// Adveccione, segun la tesis esta es la parte crucial, la que realmente hace la diferencia entre computador rectilineo y realidad curvilinea


void LatticeBoltzmann::Adveccione(void){
  int ix,iy,iz,i,rr,tt,zz;
  for(ix=0;ix<Lr;ix++)
    for(iy=0;iy<Lt;iy++)
      for(iz=0;iz<Lz;iz++)
        for(i=0;i<Q;i++){
	  rr=r(ix+V[0][i]);tt=theta(iy+V[1][i]); zz=z(iz+V[2][i]);
          if(rr<=rmax && zz<=zmax && zz>=zmin && rr>=rmin)
            f[ix+V[0][i]][iy+V[1][i]][iz+V[2][i]][i]=fnew[ix][iy][iz][i];
	  else{
	    fnew[ix][iy][iz][2]=D*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=D*f[ix][iy][iz][2];
	    fnew[ix][iy][iz][4]=D*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=D*f[ix][iy][iz][4];
	    fnew[ix][iy][iz][6]=D*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=D*f[ix][iy][iz][6];
	  }
	}
}


void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0, double Jz0){
  int ix,iy,iz,i;
  for(ix=2;ix<Lr;ix++)
    for(iy=0;iy<Lt;iy++)
      for(iz=0;iz<Lz;iz++)
        for(i=0;i<Q;i++)
          f[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i,ix);
}

// el mayor cambio viene en considerar el raiz de g que no quite de la funcion rho

void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix,iy,iz,A; double omega,rho0,Jx0,Jy0,Jz0,r0;
  omega=2*M_PI*C/lambda; A=2;
  ix=0; iy=int(thetamax/4); iz=0;
  r0=r(ix);
  rho0=A*sin(omega*t)*deltas*r0; Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i,ix);
}



/* esto tal vez lo use pero quien sabe 
void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  std::ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0,Jz0;
  for(int ix=2;ix<Lr;ix++){
    for(int iy=0;iy<Lt;iy++){
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

*/
