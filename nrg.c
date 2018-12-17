#define MAIN
#include "nrg.h"


short int u[nlinks];

int x_p[V],y_p[V],z_p[V],t_p[V],
	 x_m[V],y_m[V],z_m[V],t_m[V],v_p[V],v_m[V];
int normaener;
int offset;

int neigh_p[D],neigh_m[D];
short int stap[8];

float coseno[Nest],seno[Nest],ebco[Nest],mebco[Nest];

float resultados[NOBS_HISTER];

short int idel;

float del2,delta,cons,energia1,energia2,energiai,good;

double dener,dener2,polysumar,polysumai;
float v_dat[n_obs][maxit];

struct s_datos datos;

/*int jrand,itp,itpf,jp,jpf,indmed;*/

/*extern*/ void lee_datos(void);
/*extern*/ void lee_conf(void);
/*extern*/ void escribe_medidas(int);
/*extern*/ void escribe_conf(void);
/*extern*/ void escribe_hister(int);
/*extern*/ void tiempo(void);

/*extern*/ int Direccionamientos(void);
/*extern*/ void Inicializa(float,float,int,int);
void Inicia_hister(float,float,int,int);
/*extern*/ int Staples(int,int);
/*extern*/ void Metropolis(int,int,int);
/*extern*/ void Ajustadelta(void);
/*extern*/ void Medidas(void);

/*extern*//* int index0(int , int x[4]);*/
/*extern*//* int indexL(int, int x[4]);*/
/*extern*//* void llenax_m(int , int , int );*/

extern void Init_Rand(int);


int main(void)
{
  int ibin,j,dir,site,is,iop,nglob=L*L*L;
  int x,y,z,t,nst,it;
  float polyakov,polyr,polyi,Normaener;
  
  tiempo();
  lee_datos();

  Init_Rand(datos.seed);

  if (datos.flag>=2)
    lee_conf();
  delta=datos.delta;
  cons=(float) V*datos.mesfr;
  

  Direccionamientos();
  
  Inicializa(datos.beta,datos.gamma,datos.seed,datos.flag);
  Normaener=(float)(1./(double)normaener);

  
  for (ibin=datos.itcut;ibin<datos.nbin;ibin++)   /*   loop en numero de bloques      */
    {
      srand(datos.seed);    /* asi se reproduce la misma secuencia que
			       cuando se lee de un backup              */
      energia1=energia2=polyakov=0.0;
      
      
#ifdef HISTER
      
      
      
      /*    if(ibin<datos.nbin/2)
	{
	*/
	  datos.beta+=datos.dbeta;
	  datos.gamma+=datos.dgamma;
	  
	  /*	}
      else
	{
	  datos.beta-=datos.dbeta;
	  datos.gamma-=datos.dgamma;
	} */
      Inicia_hister(datos.beta,datos.gamma,datos.seed,datos.flag);

      
#endif
      

	for (it=0;it<datos.itmax;it++)
	  {
	    good=0.0F;
	    /*    itp=it;*/
	    for (j=0;j<datos.mesfr;j++) /*   loop de MonteCarlo sin medidas */
	      {
		/*jp=j;*/
		for (dir=0;dir<D;dir++)
		  {
		    for(site=0;site<V;site++){
		      neigh_p[4]=v_p[site];
		      neigh_m[4]=v_m[site];
		      neigh_p[3]=t_p[site];
		      neigh_m[3]=t_m[site];
		      neigh_p[2]=z_p[site];
		      neigh_m[2]=z_m[site];
		      neigh_p[1]=y_p[site];
		      neigh_m[1]=y_m[site];
		      neigh_p[0]=x_p[site];
		      neigh_m[0]=x_m[site];
		      if(neigh_p[dir])
			{
			  offset=dir*V;
			  nst=Staples(site,dir);
			  Metropolis(site,dir,nst);

			}
		    }
		  }
		
	      }
	    
	    Medidas();
	    iop=0;
	    v_dat[iop++][it] = (float) (dener *Normaener);
	    v_dat[iop++][it] = (float) (dener2 *Normaener);
	    polyr=polysumar*Normapoly;
	    polyi=polysumai*Normapoly;
	    v_dat[iop++][it] = (float) (polyr);
	    v_dat[iop++][it] = (float) (polyi);

	    
	    energia1  += dener;
	    energia2  += dener2;
	    polyakov += sqrt(polyr*polyr+polyi*polyi);
	  }
	
#ifdef HISTER
	
	resultados[0]=v_dat[0][datos.itmax-1];
	resultados[1]=v_dat[1][datos.itmax-1];
	
	
	escribe_hister(ibin);  /* fichero de histeresis */
	
#endif
	energia1 *= Normaener;
	energia2 *= Normaener;
	polyakov*= Normapoly;
	printf("\n");
	tiempo();
	printf("N=%4i, <E1>  = %+10.8f  <E2> = %+10.8f   <P>  = %10.8f",
	       ibin,energia1/datos.itmax,energia2/datos.itmax,
	       polyakov/datos.itmax);
#ifndef HISTER
	escribe_medidas(ibin);
#endif
	datos.seed=rand();
	datos.delta=delta;
	datos.itcut=ibin+1;
#ifndef HISTER
	escribe_conf();
#endif
    }
    return 1;
}

