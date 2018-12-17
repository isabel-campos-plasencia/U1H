#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mem.h>  /* No ponerla para correr en la SUN */
#include <string.h>

#define MAXIT 50
#define MAXBINS 20
#define IS_TOT 6
#define N_OBS   (2+6*IS_TOT)
#define MAXFILES 100
#define NBETAS 100

extern _stklen=10000;

float v_dat[N_OBS][MAXIT];
char *ficheros[MAXFILES];
double o0[NBETAS+1];
double o[N_OBS][NBETAS+1];

struct s_datos
{
   int itmax,		/* Numero de medidas por bloque 		 */
       mesfr,		/* frecuencia de las medidas    		 */
       nbin,		/* numero de bloque             		 */
       itcut,		/* proximo bloque a calcular    		 */
       flag,		/* conf de partida: 0(random), 1(fria),2(backup) */
       seed;		/* semilla random				 */
    float beta,		/* acoplamiento de los campos de gauge		 */
	  delta;	/* salto de Metropolis para los campos de gauge  */
};

void pinta_datos(struct s_datos *dat)
{
    printf("\n");
    printf("itmax  %d \n",dat->itmax);
    printf("mesfr  %d \n",dat->mesfr);
    printf("nbin   %d \n",dat->nbin);
    printf("itcut  %d \n",dat->itcut);
    printf("flag   %d \n",dat->flag);
    printf("seed   %d \n",dat->seed);
    printf("beta   %f \n",dat->beta);
    printf("delta  %f \n",dat->delta);
}

void lee_argumentos(int argc, char *argv[], struct s_datos *datos,
		      int *l, int *lblo, int *nblo) /*, FILE *Fnombres)*/
{
    int n1,n2,n;
    char nombres[100],name[100];
    FILE *Finput,*Fnombres;
    static struct s_datos datos_old;
    switch (argc)
    {
	case 6: sscanf(argv[6],"%s",nombres);
	case 5: sscanf(argv[4],"%d",nblo);
	case 4: sscanf(argv[3],"%d",l);
	case 3: sscanf(argv[2],"%d",&n2);
		sscanf(argv[1],"%d",&n1);  break;
	default: printf("Usage: ana n1 n2 L nblo filename\n\a");
		 exit(0);
    }
    if (n2-n1+1>MAXFILES)
    {
	printf("El numero maximo de ficheros es %i\b\n",MAXFILES);
	exit(0);
    }
    printf("\nficheros : OUT%03d.DAT - OUT%03d.DAT\n",n1,n2);
    printf("L        : %i\n",*l);

    printf("num bins : %i\n",*nblo);

    if (*nblo>MAXBINS)
    {
	printf("\nnblo no debe superar %i\n\a",MAXBINS);
	exit(0);
    }

    *lblo=(n2-n1+1)/ *nblo;
    if (*lblo==0)
	exit(0);
    n1=n2+1-(*lblo) * (*nblo);      /* para que la division en bloques sea exacta */


    if (argc==8)
    {
	Fnombres=fopen(nombres,"rt");
	if (Fnombres==NULL)
	{
	    printf("\nEl fichero %s no existe\n\a",nombres);
	    exit(0);
	}
    }
/*  else
	Fnombres=NULL;                 */

    for (n=n1;n<=n2;n++)
    {
	if (argc==8)
	    fscanf(Fnombres,"%s\n",name);
	else
	    sprintf(name,"%s%03d%s","OUT",n,".DAT");

	if( (Finput=fopen(name,"rb"))==NULL)
	{
	    printf("\nEl fichero %s no existe\n\a",name);
	    exit(0);
	}
	fread(datos,sizeof(*datos),1,Finput);
	if (n==n1)
	{
	    pinta_datos(datos);
	    memcpy(&datos_old, datos, sizeof(datos_old));
	}
	if ( (n>n1) &&
	 ( datos->beta  != datos_old.beta  ||
	   datos->itmax != datos_old.itmax ||
	   datos->mesfr != datos_old.mesfr  ))
	{
	    printf("\nLos par metros de %s no coinciden con los del anterior\n\a",name);
	    exit(0);
	}
	fclose(Finput);

	if((ficheros[n-n1]= (char *) malloc(strlen(name)))==NULL)
	{
	    printf("Fuera de memoria\b\n");
	    exit(0);
	}
	strcpy(ficheros[n-n1], name);

    }
}

void lee_datos(int n)
{
    int idat;
    FILE *Finput;
    struct s_datos dat;

    Finput=fopen(ficheros[n],"rb");
    fread(&dat,sizeof(dat),1,Finput);
    for(idat=0;idat<N_OBS;idat++)
	fread(&v_dat[idat][0],4,(size_t)dat.itmax,Finput);
    fclose(Finput);
}

void jackknife(double * p,int n, double * media, double * error)
{
    int i;
    double *sum;
    double x,Sum,sumx,sumxx;

    sum=malloc(8*n);
    memcpy(sum,p,8*n);

    Sum=0;
    for (i=0;i<n;i++)
	Sum+=sum[i];

    for (i=0;i<n;i++)
	sum[i]=Sum-sum[i];

    *media=Sum/n;

    for (i=sumx=sumxx=0;i<n;i++)
    {
	x=sum[i];
	sumx+=x;
	sumxx+=x*x;
    }
    *error=sqrt((sumxx-sumx*sumx/n)/(n*(n-1)));
}
void jackknife_der(double *pa, double *pb, double *q, int n,
			double *media, double *error)
{
    int i;
    double suma[MAXBINS],sumb[MAXBINS],sum2[MAXBINS];
    double x,Suma,Sumb,Sum2,sumx,sumxx;

/*    suma=malloc(8*n);
    sumb=malloc(8*n);
    sum2=malloc(8*n);    */


    Suma=Sumb=Sum2=0;
    for (i=0;i<n;i++)
    {
	Suma+=pa[i];
	Sumb+=pb[i];
	Sum2+= q[i];
    }
    for (i=0;i<n;i++)
    {
	suma[i]=Suma-pa[i];
	sumb[i]=Sumb-pb[i];
	sum2[i]=Sum2- q[i];
    }
    *media=Sum2/n-Suma*Sumb/(n*n);

    for (i=sumx=sumxx=0;i<n;i++)
    {
	x=sum2[i]/(n-1)-suma[i]*sumb[i]/((n-1)*(n-1));
	sumx+=x;
	sumxx+=x*x;
    }

    *error=sqrt((sumxx/n-sumx*sumx/(n*n))*(n-1)); 
/*    free(sum2);
    free(sumb);
    free(suma);*/
}

void FS(float e_min, float e_max, int numfiles,
	struct s_datos datos, float seis_v, FILE *Foutput)
{
    double delta,h,beta0,coef,expo;
    int ibeta,i,j,k;

    delta=1./(seis_v*(e_max-e_min));
    h=delta/100;
    beta0=datos.beta-delta/2.;

    for (ibeta=0;ibeta<=NBETAS;ibeta++)
    {
	o0[ibeta]=0;
	for (k=0;k<N_OBS;k++)
	    o[k][ibeta]=0;
    }

    for (i=0; i<numfiles;i++)
    {
	lee_datos(i);
	for (ibeta=0;ibeta<=NBETAS;ibeta++)
	{
	    coef=seis_v*h*(ibeta-NBETAS/2);
	    for (j=0;j<datos.itmax;j++)
	    {
		expo=exp(coef*v_dat[0][j]);
		o0[ibeta]+=expo;
		for (k=0;k<N_OBS;k++)
		    o[k][ibeta]+=v_dat[k][j]*expo;
	    }
	}
    }

    for (ibeta=0;ibeta<=NBETAS;ibeta++)
    {
	for (k=0;k<N_OBS;k++)
	    o[k][ibeta]/=o0[ibeta];
	fprintf(Foutput,"%8f %8f \n",beta0+h*ibeta,o[0][ibeta]);
    }
    for(k=20;k<=24;k=k+2)
    {
	fprintf(Foutput,"\n New plot \n");
	for (ibeta=0;ibeta<=NBETAS;ibeta++)
	{
		fprintf(Foutput,"%10e %10e \n",beta0+h*ibeta,o[k][ibeta]);
	}
		fprintf(Foutput,"\n Plot \n");
    }
}

void main(int argc, char *argv[])
{
    int l=6,
	nblo=2;

    int i,j, k,il,ib,lblo;
    float num;
    float e_min,e_max,energia;
    double sume[MAXBINS],sumee[MAXBINS];


    double v;
    double en,sigma, en_err,sigma_err;
    double corr[IS_TOT], corr_err[IS_TOT], corri[IS_TOT], corri_err[IS_TOT];
    double corrsp[IS_TOT],corrsp_err[IS_TOT];
    FILE *Foutput;

    struct s_datos datos;

    lee_argumentos(argc,argv,&datos,&l,&lblo,&nblo);
    v=pow(l,4);
    for (ib=0;ib<nblo;ib++)
    {
	sume[ib]=sumee[ib]=0;
    }
    e_min=1;
    e_max=0;
    i=-1;
    for (ib=0;ib<nblo;ib++)
    {
	for (il=0;il<lblo;il++)
	{
	    i++;
	    lee_datos(i);
/*				v_dat[0][]= Sum[Cos]/(6V)
*/
	    for (j=0;j<datos.itmax;j++)
	    {
		energia=v_dat[0][j];
		sume[ib]+=energia;
		sumee[ib]+=energia*energia;
		e_min=(energia<e_min)?energia:e_min;
		e_max=(energia>e_max)?energia:e_max;
/*        Foutput=fopen("output.dat","wa");
	fprintf(Foutput,"%i %8.6f \n",j,v_dat[12][j]); */
 
	    }
	}
    }
    num=(double) (lblo*datos.itmax);

    for (ib=0;ib<nblo;ib++)
    {
	sume[ib]/=num;
	sumee[ib]/=num;
    }
    jackknife(sume,nblo,&en,&en_err);


    printf("\n en= %8.6f(%7.6f), sigma= %8.6f(%7.6f),e î [%7.6f,%7.6f]\n",
		en,en_err,sqrt(sigma),
		sigma_err/sqrt(sigma)/2,
		e_min,e_max);

    Foutput=fopen("output.dat","wa");
}
