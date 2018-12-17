#include "nrg.h"

extern short int u[];
extern float coseno[],seno[];
extern double dener,dener2,polysumar,polysumai;
extern int neigh_p[],neigh_m[];
extern int x_p[],y_p[],z_p[],t_p[], v_p[],
		x_m[],y_m[],z_m[],t_m[],v_m[];
extern int offset;
extern short int stap[]/*,stapu[]*/;
extern int Staples(int, int);
extern int normaener;
void Medidas1(void)
{
  int dir,iblo,tirar,nst;
  int j,ja,jb,i;
  float energiaplaq,energiaplaq2;
  int x,y,z,t,is,x1,x2,x3,x4;
  int suma;

  short int uold;

  dener  = 0.0;
  dener2 = 0.0;
    
  for (dir=0;dir<D;dir++) /*   loop en direcciones de link    */
    {
      offset=dir*V;
      for(j=0;j<V;j++)
	{
	  neigh_p[4]=v_p[j];
	  neigh_m[4]=v_m[j];
	  neigh_p[3]=t_p[j];
	  neigh_m[3]=t_m[j];
	  neigh_p[2]=z_p[j];
	  neigh_m[2]=z_m[j];
	  neigh_p[1]=y_p[j];
	  neigh_m[1]=y_m[j];
	  neigh_p[0]=x_p[j];
	  neigh_m[0]=x_m[j];
	  if(neigh_p[dir])
	     {
	       nst=Staples(j,dir);

	       uold=u[j+offset];
	       for(i=0;i<nst;i++)
		 {
		   tirar=Nestm1 &
		     (stap[i]-uold);
		   energiaplaq=coseno[tirar];
		   dener+=(double) energiaplaq;
		   energiaplaq2=coseno[(2*tirar) & Nestm1];
		   dener2+=(double) energiaplaq2;
		 }
	     }
	}
    }
  
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

void Medidas0(void)
{
  int dir,nst;
  int j;



  normaener=0;
    
  for (dir=0;dir<D;dir++) /*   loop en direcciones de link    */
    {
      offset=dir*V;
      for(j=0;j<V;j++)
	{
	  neigh_p[4]=v_p[j];
	  neigh_m[4]=v_m[j];
	  neigh_p[3]=t_p[j];
	  neigh_m[3]=t_m[j];
	  neigh_p[2]=z_p[j];
	  neigh_m[2]=z_m[j];
	  neigh_p[1]=y_p[j];
	  neigh_m[1]=y_m[j];
	  neigh_p[0]=x_p[j];
	  neigh_m[0]=x_m[j];
	  if(neigh_p[dir])
	     {
	       nst=Staples(j,dir);
	       normaener+=nst;

	     }
	}
    }
  /*  printf("\n normaener = %d",normaener);*/
}




