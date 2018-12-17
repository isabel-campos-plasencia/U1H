#define MAIN
#include "nrg.h"

short int u[nlinks];

int x_p[L],y_p[L],z_p[L],t_p[L],
    x_m[L],y_m[L],z_m[L],t_m[L];



int neigh_p[4],neigh_m[4];
int stap[6];

float coseno[Nest],seno[Nest],ebco[Nest],mebco[Nest];

float resultados[NOBS_HISTER];

int idel,offset;

float del2,delta,cons,energia1,energia2,energiai,good,tot,aceptancia;

double dener,dener2,polysumar,polysumai;
float v_dat[n_obs][maxit];

struct s_datos datos;

extern void lee_datos(void);
extern void lee_conf(void);
extern void escribe_medidas(int);
extern void escribe_conf(void);
extern void escribe_hister(int);
extern void tiempo(void);

extern void Direccionamientos(void);
extern void Inicializa(float,float,int,int);
extern void Inicia_hister(float,float,int,int);
extern void Staples(int,int);
extern void Metropolis(int,int);
extern void Ajustadelta(void);
extern void Medidas(void);
extern void Init_Rand(int);


int main(void)
{
	int ibin,it,j,dir,site,is,iop,nglob=L*L*L;
    int x,y,z,t;
    float polyakov,polyr,polyi;

    tiempo();
    lee_datos();


    Init_Rand(datos.seed);

    if (datos.flag>=2)
	lee_conf();
    delta=datos.delta;
	cons=(float) V*datos.mesfr;

    Direccionamientos();
    Inicializa(datos.beta,datos.gamma,datos.seed,datos.flag);


    for (ibin=datos.itcut;ibin<datos.nbin;ibin++)   /*   loop en numero de bloques      */
    {
	srand(datos.seed);    /* asi se reproduce la misma secuencia que
				 cuando se lee de un backup              */
	energia1=energia2=0.0;


#ifdef HISTER

	

	if(ibin<datos.nbin/2)
	  {     
	    datos.beta+=datos.dbeta;
	    datos.gamma+=datos.dgamma;
	    
	  }  
	else
	  {
	    datos.beta-=datos.dbeta;
	    datos.gamma-=datos.dgamma;
	  }  
        Inicia_hister(datos.beta,datos.gamma,datos.seed,datos.flag);


#endif
	polyakov=0.0;

	for (it=0;it<datos.itmax;it++)
	{
	    tot=good=0.0F;
	    for (j=0;j<datos.mesfr;j++) /*   loop de MonteCarlo sin medidas */
	    {

		for (dir=0;dir<4;dir++)
		{
		    site=0;
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
				    Staples(site,dir);
				    Metropolis(site,dir);
				    site++;
				}
			    }
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
	aceptancia = 100.*good/tot;

	energia1 *= Normaener;
	energia2 *= Normaener;
	polyakov*= Normapoly;
	printf("\n");
	tiempo();
	printf("N=%4i, Acep = %f  <E1>  = %+10.8f  <E2> = %+10.8f   <P>  = %10.8f",
			ibin,aceptancia,energia1/datos.itmax,energia2/datos.itmax,
	       polyakov/datos.itmax);
#ifndef HISTER
	escribe_medidas(ibin);
#endif                    
	datos.seed=RANDOM;
	datos.delta=delta;
	datos.itcut=ibin+1;
#ifndef HISTER
	escribe_conf();
#endif

	}
    return 1;
}
