#include "nrg.h"


/*extern int jrand;*/

extern short int u[];
extern float coseno[],seno[];
extern double dener,dener2,polysumar,polysumai;
/*extern int neigh_p[],neigh_m[];*/
extern int x_p[],y_p[],z_p[],t_p[], v_p[],
		x_m[],y_m[],z_m[],t_m[],v_m[];

extern short int stap[]/*,stapu[]*/;


void Medidas(void)
{
  int dir,offset,ortog,offortog,iblo,tirar;
  int j,ja,jb,neigb_p[D],neigh_p[D];
  float energiaplaq,energiaplaq2;
  int x,y,z,t,is,x1,x2,x3,x4;
  int suma;
  dener  = 0.0;
  dener2 = 0.0;
  
  for(j=0;j<V;j++){
    neigh_p[4]=v_p[j];
    neigh_p[3]=t_p[j];
    neigh_p[2]=z_p[j];
    neigh_p[1]=y_p[j];
    neigh_p[0]=x_p[j];
    for (dir=0;dir<D-1;dir++) /*   loop en direcciones de link    */
      {
	offset=dir*V;
	if(neigh_p[dir]){
	  jb=j+neigh_p[dir];
	  neigb_p[4]=v_p[jb];
	  neigb_p[3]=t_p[jb];
	  neigb_p[2]=z_p[jb];
	  neigb_p[1]=y_p[jb];
	  neigb_p[0]=x_p[jb];
	  for (ortog=dir+1;ortog<D;ortog++)
	    {
	      if(neigh_p[ortog]&&neigb_p[ortog])
		{
		  offortog=ortog*V;
		  ja=j+neigh_p[ortog];
		  tirar=Nestm1 &
		    (u[j+offset]+u[jb+offortog]
		     -u[ja+offset]-u[j+offortog]);
		  energiaplaq=coseno[tirar];
		  /*
if(j==jrand)
  {
printf("\n Medidas \n j=%d, dir=%d, offset = %d, u[j,dir] =%d",j,dir,offset,u[j+offset]);
printf("\n ortog=%d,ja=%d,jb=%d,offortog = %d,-u[j+ortog]=%d,-u[ja]=%d,u[jb]=%d",ortog,ja,jb,offortog,-u[j+offortog],-u[ja+offset],u[jb+offortog]);
printf("\n thetapalaq = %d, energiaplaq = %f",tirar,energiaplaq);
  }*/
		  dener+=(double) energiaplaq;
		  energiaplaq2=coseno[(2*tirar) & Nestm1];
		  dener2+=(double) energiaplaq2;
		}
	    }  /*for ortog*/
	}   /*if neigh_p[dir]*/
      }  /*for dir*/
  }  /*for j*/
  /*  if(j==jrand)printf("\n energia = %f",dener);  */
  polysumar=0.0;
  polysumai=0.0;
  
  /*  for (dir=0;dir<D;dir++)
      for(is=0;is<V;is++){
          for (x1=0;x1<L;x1++)
	    for (x2=0;x2<L;x2++)
	    for (x3=0;x3<L;x3++)
	    {
	    suma=0;
	    for (x4=0;x4<L;x4++)
	    {
	    switch(dir)
	    {
	    case 0:
	    is=x4+x1*L+x2*L*L+x3*L*L*L;
	    break;
	    case 1:
	    is=x1+x4*L+x2*L*L+x3*L*L*L;
	    break;
	    case 2:
	    is=x1+x2*L+x4*L*L+x3*L*L*L;
	    break;
	    case 3:
	    is=x1+x2*L+x3*L*L+x4*L*L*L;
	    break;
	    }       
      suma+=u[is+dir*V];
      polysumar+=coseno[suma&Nestm1];
      polysumai+=seno[suma&Nestm1];
    }*/
  
}





