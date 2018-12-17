#include "nrg.h"

extern int x_p[],y_p[],z_p[],t_p[],
	   x_m[],y_m[],z_m[],t_m[];


char dir[60];
extern struct s_datos datos;


extern short int u[];

extern float coseno[],seno[],ebco[],mebco[];
extern float del2,delta,resultados[];
extern int idel;
extern int offset;



void Init_Rand(int semilla)     /* Se utiliza el generador estandar de C  */
{                               /* para inicializar otro generador random */
    int i;                      /* mas seguro. Para aumentar mas la segu- */
                                /* ridad se eliminan los primeros numeros */
    srand((unsigned)semilla);
    for (i=0;i<111;i++)         /* Probablemente no es necesario */
        rand();

    ip=128;    
    ip1=ip-24;    
    ip2=ip-55;    
    ip3=ip-61;
    
    for (i=0; i<256; i++)
        ira[i] = (unsigned) rand()+ (unsigned) rand(); 

    for (i=0;i<1111;i++)        /* Probablemente tampoco es necesario */
        RANDOM;
}



void Direccionamientos(void)
{
    int i;

    for (i=0;i<L;i++)
    {
	x_p[i]= 1;
	x_m[i]=-1;
	y_p[i]= L;
	y_m[i]=-L;
	z_p[i]= L*L;
	z_m[i]=-L*L;
	t_p[i]= L*L*L;
	t_m[i]=-L*L*L;
    }
    x_m[0]= L-1;
    y_m[0]=(L-1)*L;
    z_m[0]=(L-1)*L*L;
    t_m[0]=(L-1)*L*L*L;
    x_p[L-1]=-x_m[0];
    y_p[L-1]=-y_m[0];
    z_p[L-1]=-z_m[0];
    t_p[L-1]=-t_m[0];
}

void Inicializa(float beta, float gamma, int semilla, int flag)
{
    int i,x,y,z,t,dir;
    int igen;
    double theta,dcos,debco;

    srand(semilla);

    for (i=0;i<Nest;i++)
    {
	theta = (double) i * twopioverNest;
	dcos  = cos(theta);
	debco = exp(beta*dcos + gamma*cos(2.0*theta));
	ebco[i]   = (float) debco;
	mebco[i]  = (float) (1.0/debco);
	coseno[i] = (float) dcos;
	seno[i]   = (float) sin(theta);
    }
    idel=(int) (delta/twopioverNest);
    del2= 2.0F * (float) idel + 1.0F;

    if (flag!=2)
    {
	igen=0;
	for (dir=0;dir<4;dir++)
	    for (t=0;t<L;t++)
		for (z=0;z<L;z++)
		    for (y=0;y<L;y++)
			for (x=0;x<L;x++)
			{

			    if (flag == 0)
			      {
				u[igen++]= (int) (FRANDOM*rNest);
			      }
			    else
				u[igen++]=0;
			}
    }


}

void Inicia_hister(float beta, float gamma, int semilla, int flag)
{
    int i,x,y,z,t,dir;
    int igen;
    double theta,dcos,debco;

    srand(semilla);

    for (i=0;i<Nest;i++)
    {
	theta = (double) i * twopioverNest;
	dcos  = cos(theta);
	debco = exp(beta*dcos + gamma*cos(2.0*theta));
	ebco[i]   = (float) debco;
	mebco[i]  = (float) (1.0/debco);
	coseno[i] = (float) dcos;
	seno[i]   = (float) sin(theta);
    }
    idel=(int) (delta/twopioverNest);
    del2= 2.0F * (float) idel + 1.0F;


}






void escribe_hister(int i)
{

    FILE *Fhister;
    char name[100];
    int k;


    sprintf(name,"%s%s_%d",dir,"hister",L);

    Fhister=fopen(name,"a");    
    fprintf(Fhister,"%7.4f %7.4f ",datos.beta,datos.gamma);
    for (k=0;k<NOBS_HISTER;k++)
        fprintf(Fhister,"%f ",resultados[k]);
    fprintf(Fhister,"\n");    
    fclose(Fhister);   

    printf("(%5.3f,%5.3f),",datos.beta,datos.gamma);

}





