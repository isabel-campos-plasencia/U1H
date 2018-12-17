#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


#define MAXIT 50
#define MAXBINS 20

#define N_OBS   2
#define MAXFILES 1000
#define NBETAS 20
#define NDISCR 100



float v_dat[N_OBS][MAXIT];
char *ficheros[MAXFILES];
double o0[NBETAS+1];
double o[N_OBS][NBETAS+1];
double oE[N_OBS][NBETAS+1];
int frec[MAXBINS][NDISCR+1];

struct s_datos
{
int itmax,		/* Numero de medidas por bloque 		 */
       mesfr,		/* frecuencia de las medidas    		 */
       nbin,		/* numero de bloque             		 */
       itcut,		/* proximo bloque a calcular    		 */
       flag,		/* conf de partida: 0(random), 1(fria),2(backup) */
       seed;		/* semilla random				 */
    float beta,		/* acoplamiento de los campos de gauge		 */
          gamma,        /* acoplamiento de plaquetas cuadrado            */
	  delta;	/* salto de Metropolis para los campos de gauge  */
};
void Histog(float min, float max, int nblo,FILE *Foutput)
{
int suma,ie,ib;
float energ;
fprintf(Foutput,"new plot \n set limits x %f %f \n",min,max);
	for (ie=0;ie<NDISCR;ie++)
	{
	suma=0;
	energ=((float)ie*(max-min)/(float)NDISCR)+min;
	for (ib=0;ib<nblo;ib++)
		suma+=frec[ib][ie];
	fprintf(Foutput,"%f %d\n",energ,suma);
	}
	fprintf(Foutput,"bargraph\n");
}




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
    printf("gamma   %f \n",dat->gamma);
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
	case 6: sscanf(argv[5],"%s",nombres);
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


    if (argc==6)
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
	if (argc==6)
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
	   datos->gamma  != datos_old.gamma  ||
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
	x=sum2[i]/(n-1)-suma[i]*sumb[i]/(n-1)/(n-1);
	sumx+=x;
	sumxx+=x*x;
    }

    *error=sqrt((sumxx/n-sumx*sumx/(n*n))*(n-1));
/*    free(sum2);
    free(sumb);
    free(suma);*/
}

void FS(float sigma, int numfiles,
	struct s_datos datos, float seis_v, FILE *Foutput)
{
    double delta,h,beta0,coef,expo;
    int ibeta,i,j,k;

    delta=2./(seis_v*sigma);
    h=delta/NBETAS;
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
		{
		    o[k][ibeta]+=v_dat[k][j]*expo;
            oE[k][ibeta]+=v_dat[k][j]*v_dat[0][j]*expo;
		}
	    }
	}
    }

    for (ibeta=0;ibeta<=NBETAS;ibeta++)
    {
	for (k=0;k<N_OBS;k++)
	{
	    o[k][ibeta]/=o0[ibeta];
	    oE[k][ibeta]/=o0[ibeta];
	}
	fprintf(Foutput,"%8f %8f ",beta0+h*ibeta,o[0][ibeta]);
	fprintf(Foutput,"\n");
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
    double en,sigma2, en_err,sigma2_err;
	FILE *Foutput,*fevol;
	float h,c1,c2,sum_par;
	int ie,ii;

    struct s_datos datos;

#ifdef GRAF
    int g_driver=VGA,g_mode=2;
    float red;
#endif

    lee_argumentos(argc,argv,&datos,&l,&lblo,&nblo);
    v=pow(l,4);
    for (ib=0;ib<nblo;ib++)
    {
	sume[ib]=sumee[ib]=0;
	}
    for(ib=0;ib<nblo;ib++)
	for (i=0;i<=NDISCR;i++)
	    frec[ib][i]=0;
    e_min=1;
    e_max=0;
	i=-1;
	fevol=fopen("mcevol.plt","w");
	for (ib=0;ib<nblo;ib++)
	{
	for (il=0;il<lblo;il++)
	{
	    i++;
	    lee_datos(i);
	    for (j=0;j<datos.itmax;j++)
	    {
		energia=v_dat[0][j];
		e_min=(energia<e_min)?energia:e_min;
		e_max=(energia>e_max)?energia:e_max;
	    }
	}
    }

    h=1.0/(e_max-e_min);
	fprintf(fevol,"set limits y %f %f \n set limits x 0 %d \n",e_min,e_max,
				   nblo*lblo*datos.itmax);
#ifdef GRAF
	initgraph(&g_driver, &g_mode, "c:\\borlandc\\bgi");
	red=10.;
#endif

	i=-1;
	ii=0;
	for (ib=0;ib<nblo;ib++)
	{
	for (il=0;il<lblo;il++)
	{
	    i++;
		lee_datos(i);

/*				v_dat[0][]= Sum[Cos]/(6V) */
		for (j=0,sum_par=0.;j<datos.itmax;j++)
		{
		energia=v_dat[0][j];
		sum_par+=energia;
		sume[ib]+=energia;
		sumee[ib]+=energia*energia;
		ie=(int) ((v_dat[0][j]-e_min)*h*NDISCR);
		frec[ib][ie]++;
		fprintf(fevol," %d %f \n",ii++,energia);
		if(ii%512==0)fprintf(fevol,"plot\n");
		}
	}
    fprintf(fevol," plot \n");
#ifdef GRAF
	plotfrec(&red,ib+1);
#endif
    }
#ifdef GRAF
    printf("\a");
    getche();
    closegraph();
#endif

    num=(double) (lblo*datos.itmax);

    for (ib=0;ib<nblo;ib++)
    {
	sume[ib]/=num;
	sumee[ib]/=num;

    }
    jackknife(sume,nblo,&en,&en_err);
    jackknife_der(sume,sume,sumee,nblo,&sigma2,&sigma2_err);


    printf("\n en= %8.6f(%7.6f), sigma= %8.6f(%7.6f),e î [%7.6f,%7.6f]\n",
		en,en_err,sqrt(sigma2),
		sigma2_err/sqrt(sigma2)/2,
		e_min,e_max);



	fclose(fevol);
	Foutput=fopen("histog.dat","w");
	Histog(e_min,e_max,nblo,Foutput); 
    fclose(Foutput);
/*	Foutput=fopen("output.dat","wt");
	FS(sqrt(sigma2),nblo*lblo,datos,6*v,Foutput);
	fclose(Foutput);
*/
}

