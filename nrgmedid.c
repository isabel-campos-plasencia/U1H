#include "nrg.h"

extern short int u[];
extern float coseno[],seno[];
extern double dener,dener2,polysumar,polysumai;
extern int neigh_p[],neigh_m[];
extern int x_p[],y_p[],z_p[],t_p[],
	   x_m[],y_m[],z_m[],t_m[];

extern int stap[],stapu[];


void Medidas(void)
{
    int dir,offset,ortog,offortog,iblo,tirar;
    int j,ja,jb;
    float energiaplaq,energiaplaq2;
    int x,y,z,t,is,x1,x2,x3,x4;
    int suma;
    dener  = 0.0;
    dener2 = 0.0;

   for (dir=0;dir<4;dir++) /*   loop en direcciones de link    */
    {
	j=0;
	offset=dir*V;
	for (t=0;t<L;t++)
	{
	    neigh_p[3]=t_p[t];
	    neigh_m[3]=t_m[t];
	    for (z=0;z<L;z++)
	    {
		neigh_p[2]=z_p[z];
		neigh_m[2]=z_m[z];
		for (y=0;y<L;y++)
		{
		    neigh_p[1]=y_p[y];
		    neigh_m[1]=y_m[y];
		    for (x=0;x<L;x++)
		    {
			neigh_p[0]=x_p[x];
			neigh_m[0]=x_m[x];
			jb=j+neigh_p[dir];
			for (ortog=0;ortog<4;ortog++)
			{
			    if (ortog>dir)
			    {
				offortog=ortog*V;
				ja=j+neigh_p[ortog];
				tirar=Nestm1 &
					      (u[j+offset]+u[jb+offortog]
					      -u[ja+offset]-u[j+offortog]);

				energiaplaq=coseno[tirar];
				dener+=(double) energiaplaq;
				energiaplaq2=coseno[(2*tirar) & Nestm1];
				dener2+=(double) energiaplaq2;
			    }
			}
			j++;
		    }
		}
	    }
	}
    }

  polysumar=0.0;
  polysumai=0.0;

  for (dir=0;dir<4;dir++)
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
	    }
	  polysumar+=coseno[suma&Nestm1];
	  polysumai+=seno[suma&Nestm1];
	  }

}

