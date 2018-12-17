#include "nrg.h"

char dir[60];
extern short int u[];
extern float v_dat[6][maxit];
extern struct s_datos datos;

FILE *Finput,*Foutput,*Fconfig;


void tiempo(void)
{
    static time_t time1=0,
                  time2;
    static clock_t clock1=0,
                   clock2;
    int temp;
    float t_CPU;    

    time2 = time(NULL);
    clock2 = clock();
    if (time1>0)
    { 
        temp=time2-time1;
        t_CPU=(0.5+(float)(clock2-clock1))/CLOCKS_PER_SEC;        
        if (temp)
            printf("%4ds (%2.0f%%): ",temp,t_CPU/(temp+0.5)*100.);
        else
            printf("%4ds (??%%): ",temp);
    }
    time1=time2;
    clock1=clock2;
    
}
void pinta_datos(struct s_datos *dat)   /* solo para debug */
{
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
void lee_datos(void)
{
    int j;
    int * ptdatos_int;
    float * ptdatos_real;
    char name[60];


    Finput=fopen("input","r");
    if (Finput==NULL)
    {
	printf(" No existe el fichero INPUT.\n");
	exit(0);
    }
    fscanf(Finput,"%s",dir);
    for (j=0,ptdatos_int=&datos.itmax;j<NDATOS_INT;j++)
	fscanf(Finput,"%d",ptdatos_int++);
    for (j=0,ptdatos_real=&datos.beta;j<NDATOS_FLOAT;j++)
	fscanf(Finput,"%f",ptdatos_real++);
    fclose(Finput);
#ifdef DEBUG
    pinta_datos(&datos);
#endif
    if (datos.itmax>maxit)
    {
	printf(" itmax > %i\a\n",maxit);
	exit(0);
    }
    if (datos.flag<2)
    {                                      /* existe outxxx.dat o conf. ? */

	sprintf(name,"%s%s%03d%s",dir,"OUT",datos.itcut,".DAT");
	Foutput=fopen(name,"rb");
	if (Foutput!=NULL)
	{
	    fclose(Foutput);
	    printf(" %s  ya existe.\a\n",name);
	    /* exit(0); */
	}
	sprintf(name,"%s%s",dir,"conf");

	Foutput=fopen(name,"rb");
	if (Foutput!=NULL)
	{
	    fclose(Foutput);
	    printf(" %s  ya existe.\a\n","CONF.");
       /*	    exit(0);        */
	}
    }
}

void lee_conf(void)
{
    struct s_datos datosb;
    char name[60];

    sprintf(name,"%s%s",dir,"conf");

    Fconfig=fopen(name,"rb");
    if (Fconfig==NULL)
    {
	printf(" No existe el fichero CONF.\a\n");
	exit(0);
    }

    fread(&datosb,sizeof(datosb),1,Fconfig);
    fread(u,sizeof(short int),nlinks,Fconfig);
#ifdef DEBUG
    pinta_datos(&datosb);
    printf("%4d %4d %4d\n",interior[0],interior[1],interior[nlinks-1]);
#endif
    if (datos.itmax != datosb.itmax ||
	datos.mesfr != datosb.mesfr ||
	datos.beta  != datosb.beta  ||
	datos.gamma != datosb.gamma ||
	datos.flag  == 3 )
    {
	printf(" Los datos de CONF. no son compatibles con los de INPUT.\n");
	printf(" Se utilizaran solo los de INPUT.\a\n");

	sprintf(name,"%s%s%03d%s",dir,"OUT",datos.itcut,".DAT");

	Foutput=fopen(name,"rb");
	if (Foutput!=NULL)
	{
	    fclose(Foutput);
	    printf(" %s  ya existe.\a\n",name);
	    exit(0);
	}
    }
    else
    {
         datos.itcut=datosb.itcut; 
	datos.delta=datosb.delta;
	datos.seed=datosb.seed;
    }
    fclose(Fconfig);
}


void escribe_medidas(int i)
{
    int idat;
    char name[60];

    sprintf(name,"%s%s%03d%s",dir,"OUT",i,".DAT");

    Foutput=fopen(name,"wb");
    fwrite(&datos,sizeof(datos),1,Foutput);
    for(idat=0;idat<n_obs;idat++)
	fwrite(&v_dat[idat][0],4*datos.itmax,1,Foutput);
    fclose(Foutput);
}

void escribe_conf(void)
{
    char name[60],name_dollar[60],name_old[60];
    sprintf(name_dollar,"%s%s",dir,"conf.$$$");
    sprintf(name,"%s%s",dir,"conf");
    sprintf(name_old,"%s%s",dir,"conf.old");

    Fconfig=fopen(name_dollar,"wb");
    fwrite(&datos,sizeof(datos),1,Fconfig);
#ifdef DEBUG
    pinta_datos(&datos);
#endif
    fwrite(u,sizeof(short int),nlinks,Fconfig);

    fclose(Fconfig);
    remove(name_old);
    rename(name,name_old);
    rename(name_dollar,name);
}


