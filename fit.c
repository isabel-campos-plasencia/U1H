/*
			Programa de analisis ana_tot.c  V2.0  (5/12/93)


Analiza ficheros de datos procedentes de los programas u1h, uqh, nrg, lam, etc.

Produce los siguientes graficos:

		-  Evolucion de la energia (I columna de resultados) en #MC
		-  Histograma de la energia y error
		-  Energia + observable(idat) en funcion de beta
		-  Derivadas de los anteriores
		-  Ajuste a cubicas de los 2 picos del histograma


		GRAF: dibuja histograma durante el analisis 
                      (solo para BORLANDC o gcc bajo X11)


Comentarios:



*/

#define MANUAL    /*ajuste manual de la anchura de los picos*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


#define maxit 500
#define nbetas 200
#define nener 100
#define maxbins 20
#define numdat 8*4     /* maximo 8 Observables x 4 RGT */ 
#define NP 4
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}



double a[NP][NP],b[NP];

float v_dat[numdat][maxit];
long int frec[maxbins][nener+1];
double xfrec[maxbins][nener+1];
double coup_maxder[maxbins+1],maxder[maxbins+1];
double en_v[nbetas+1],der_v[nbetas+1],en_err[nbetas+1],der_err[nbetas+1];
double e0,coupn,shift,e[nener+1],f[nener+1];
long numit,ittot,itfrec;
double O_v[nbetas+1],O_err[nbetas+1],derO_v[nbetas+1],derO_err[nbetas+1];

struct s_datos
{
long int itmax,		/* Numero de medidas por bloque 		 */
       mesfr,		/* frecuencia de las medidas    		 */
       nbin,		/* numero de bloque             		 */
       itcut,		/* proximo bloque a calcular    		 */
       flag,		/* conf de partida: 0(random), 1(fria),2(backup) */
       seed;		/* semilla random				 */
    float kappa1,		/* coupling 1os vecinos */
	  kappa2,		/* coupling 2os vecinos */
          beta;                 /* coupling gauge       */ 
};

char formato[]="new plot\n"
	       "stamp top date\n"
	       "%s\n%s\n"
	       "title left font T '%s'\n"
	       "case              '%s'\n"
	       "title bottom font T '%s'\n"
	       "case                '%s'\n"
	       "set order x y dy\n"
	       "set labels left font T\n"
	       "set labels bottom font T\n"
	       "set symbol %s\n";

void gaussj(void)
{
    int indxc[NP],indxr[NP],ipiv[NP];
    int i,icol,irow,j,k,l,ll;
    double big,dum,pivinv;
    int n=NP;

    for (j=0;j<n;j++) ipiv[j]=0;

    for (i=0;i<n;i++) 
    {
	big=0.0;
	for (j=0;j<n;j++)
        {
	    if (ipiv[j] != 1)
		for (k=0;k<n;k++)
		{
		    if (ipiv[k] == 0)
		    {
			if (fabs(a[j][k]) >= big) 
			{
			    big=fabs(a[j][k]);
			    irow=j;
			    icol=k;
			}
		    } 
		    else if (ipiv[k] > 1)
		    {
			printf("\a Matriz singular\n");
			exit(0);
		    }
		}
        }
	++(ipiv[icol]);
	if (irow != icol) 
	{
	    for (l=0;l<n;l++)
            {
		SWAP(a[irow][l],a[icol][l])
		SWAP(b[irow],b[icol])
            }
	}
	indxr[i]=irow;
	indxc[i]=icol;
	if (a[icol][icol] == 0.0)
	{
	    printf("\a Matriz singular\n");
	    exit(0);
	}
	pivinv=1.0/a[icol][icol];
	a[icol][icol]=1.0;
	for (l=0;l<n;l++) a[icol][l] *= pivinv;
	    b[icol] *= pivinv;
	for (ll=0;ll<n;ll++)
	    if (ll != icol)
	    {
		dum=a[ll][icol];
		a[ll][icol]=0.0;
		for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
			b[ll] -= b[icol]*dum;
	    }
    }
    for (l=n-1;l>=0;l--) 
    {
        if (indxr[l] != indxc[l])
	    for (k=1;k<=n;k++)
	        SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }

}

void cubic_fit(double *e,double *f,
		double frac, int ndata,double *vsol,double *vb, int *used)
{
    double maximo,ancho;
    int i,j,k,u1,u2,u3,u4;
    double discri,d2,sol1,sol2,pj;
    int imax[2];
    double sum[7],sumf[4];
    if (frac) 
    {
	if (vsol[0]==0)
	{
	    maximo=0;
	    for (i=0;i<ndata/2;i++)
		if (f[i]>maximo)
		{
		    maximo=f[i];
		    imax[0]=i;
		}
	    maximo=0;
	    for (i=ndata-1;i>ndata/2;i--)
		if (f[i]>maximo)
		{
		    maximo=f[i];
		    imax[1]=i;
		}
	    for(i=0;i<2;i++)
		vsol[i]=e[imax[i]];
	}
	ancho=(vsol[1]-vsol[0])/frac;
#ifndef MANUAL
	for (i=0;i<ndata;i++)
	{
	    if (e[i]<vsol[0]-ancho) used[0]=i+1;
	    if (e[i]<vsol[0]+ancho) used[1]=i;
	    if (e[i]<vsol[1]-ancho) used[2]=i+1;
	    if (e[i]<vsol[1]+ancho) used[3]=i;
	}
#endif    
    }
    for (k=0;k<2;k++)
    {
	for (i=0;i<7;i++)
	    sum[i]=0;
        for (i=0;i<4;i++)
	    sumf[i]=0;

	for (i=used[2*k];i<=used[2*k+1];i++)
	{
            pj=1;
	    for (j=0;j<7;j++)
            {
		sum[j] += pj;
                pj*=e[i];
	    }
        }

	for (i=used[2*k];i<=used[2*k+1];i++)
	{
            pj=f[i];
	    for (j=0;j<4;j++)
            {
		sumf[j] += pj;
		pj*=e[i];
	    }
        }
        for (i=0;i<4;i++)
	    for(j=0;j<4;j++)
            {
		a[i][j]=sum[i+j];
	        b[i]=sumf[i];
	    }
	gaussj();

	d2=b[2]*b[2]-3*b[3]*b[1];
	if (d2<=0)
	{
	     printf("El discriminante es %lf \n",d2);
	     exit(0);
	}
	discri=sqrt(d2);
	sol1=(-b[2]-discri)/b[3]/3;
	sol2=(-b[2]+discri)/b[3]/3;
	if (2*b[2]+6*b[3]*sol1>0)
	{
	    SWAP(sol1,sol2)
	}
	vsol[k]=sol1;
	for (i=0;i<4;i++)
	    vb[4*k+i]=b[i];
    }
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
    printf("kappa1   %f \n",dat->kappa1);
    printf("kappa2   %f \n",dat->kappa2);
    printf("beta  %f \n",dat->beta); 
}


void titulo(char *nombre, char *nombre_case, struct s_datos *dat,
	      int n1, int n2, int l, int nblo)
{
    int i;


    sprintf(nombre,"title top font T size 3 'L=%i, "
		    "#MC=%lix%lix(%i-%i), NBbE=%i'",
		    l,dat->itmax,dat->mesfr,n1,n2,nblo);

    memcpy(nombre_case,"case                    ",24);
    i=23;
    do
    {
	switch (nombre[++i])
	{
	    case 39 :
	    case 0  : nombre_case[i]=nombre[i]; break;
	    case 'b': if (nombre[i+1]!='=') { nombre_case[i]=' ';break;}
	    case 'k': nombre_case[i]='G'; break;
	    case 'l': nombre_case[i]='G';break;
            case 'x': nombre_case[i]='S'; break;
	    case 'B':
	    case 'E': nombre_case[i]='-'; break;
	    default: nombre_case[i]=' ';
	}
    }
    while (nombre[i]!=0);
}


void main(int argc, char *argv[])
{
    char nombre[100],nombre_case[100];
    int i,j,k,n,n1,n2,idat,jdat_c;
    int lblo,ib,ibb,it,jdat;
    long v;
    double datos_coupling[3];
    int vol_coupling[3];

    int l   =6,       /* default  */
	nblo=2;
    /* variables del cubic-fit*/
    double frac=0;
    double e0,coupn,shift,vsol[2],e[nener+1],f[nener+1],xx;
    double vb[8],fr[2];
    int ndata,used[4];
    double E1,E2,sumC,sumCC,sumE[2],sumEE[2],E_err[2],coup_b,coupn_err;
    double suma,cte,yy,FR1,FR2,sumL,sumLL,lat_err;

    char name[12],n_coup[12],n_ener[12],n_deri[12];
    FILE *Finput,*Foutplt,*Fnombres,*Fh;
    struct s_datos datos,datos_old;

    double ymin,ymax,expo_frec;
    double coup,c1,c2,delta;
    double h,x,y,x0,y0,en,sigma,sum,sume,sume2,cv;
    double Sum,Sum2,Sumder,Sumder2;
    double derivada,err,err_maxder;
    long int sumf;

    double sumf2,media,disp,hist[nener+1],normalizacion;

    char nombres[100];


    double SumO,SumO2,SumderO,SumderO2;
    double sumO,sumOe,sumx;
    double expo,O,derO;

    switch (argc)        /* n1 n2 l nblo jdat frac nombres */
    {
    case 8: sscanf(argv[7],"%s",nombres);
    case 7: sscanf(argv[6],"%lf",&frac);
    case 6: sscanf(argv[5],"%d",&jdat);
    case 5: sscanf(argv[4],"%d",&nblo);
    case 4: sscanf(argv[3],"%d",&l);
    case 3: sscanf(argv[2],"%d",&n2);
      sscanf(argv[1],"%d",&n1);  break;
    default: printf("Usage: ana_fs n1 n2 L nblo idat frac filename\n\a");
      exit(0);
    }
    printf("\nficheros : OUT%03d - OUT%03d\n",n1,n2);
    printf("L        : %i\n",l);

    printf("num bins : %i\n",nblo);

    if (nblo>maxbins)
    {
	printf("\nnblo no debe superar %i\n\a",maxbins);
	exit(0);
    }


    lblo=(n2-n1+1)/nblo;
    if (lblo==0)
	exit(0);
    n1=n2+1-lblo*nblo;       /* para que la division en bloques sea exacta */

    memset(&datos_old,0,sizeof(datos_old));
    sum=sume=sume2=0.0;
    ymin=10;
    ymax=-10;
			       /* hace una primera lectura para calcular */
			       /* los limites, dispersion y media        */

    if (argc==8)
    {
        Fnombres=fopen(nombres,"rt");
        if (Fnombres==NULL)
        {
            printf("\nEl fichero %s no existe\n\a",nombres);
            exit(0);
        }
    }



    for (n=n1;n<=n2;n++)
    {

      if (argc==8)
	fscanf(Fnombres,"%s\n",name);
      else
	sprintf(name,"%s%03d%s","OUT",n,".DAT");

      Finput=fopen(name,"rb");
      if (Finput==NULL)
	{
	  printf("\nEl fichero %s no existe\n\a",name);
	  exit(0);
	}
      fread(&datos,sizeof(datos),1,Finput);
      
  
      if (n==n1) pinta_datos(&datos);
      
      if ( (n>n1) &&
	   (datos.kappa1  != datos_old.kappa1  ||
	    datos.kappa2  != datos_old.kappa2  ||
	    datos.beta  != datos_old.beta  ||
	    datos.itmax != datos_old.itmax ||
	    datos.mesfr != datos_old.mesfr  ))
	{
	  printf("\nLos parametros de %s son diferentes\n\a",name);
	  exit(0);
	}
      memcpy(&datos_old, &datos, sizeof(datos));
      
      for(idat=0;idat<jdat+1;idat++)
	fread(&v_dat[idat][0],4,(size_t)datos.itmax,Finput);

      fclose(Finput);
      
      for (it=0;it<datos.itmax;it++)
	{
	  y=v_dat[jdat][it];
	  sum++;
	  sume += y ;
	  sume2+= y*y;
	  ymin= (y<ymin)?y:ymin;
	  ymax= (y>ymax)?y:ymax;
	}
    }
    if (argc==8)
      fclose(Fnombres);
    v=(long)l*l*l*l;

    vol_coupling[0] = 6*v;
    vol_coupling[1] = 4*v;
    vol_coupling[2] = 6*v;

    datos_coupling[0] = datos.beta;
    datos_coupling[1] = datos.kappa1;
    datos_coupling[2] = datos.kappa2;

    v=vol_coupling[jdat];

    coup=datos_coupling[jdat];


    strcpy(n_ener,"EBplaqE");
    strcpy(n_deri,"EBplaqE,b");
    strcpy(n_coup,"b");

    en=sume/sum;
    cv=sume2/sum-en*en;
    sigma=sqrt(cv);
    printf("\nen= %10f, sigma= %10f, en_min = %10f, en_max= %10f\n\n",
		en,sigma,ymin,ymax);

    for(ib=0;ib<nblo;ib++)
	for (i=0;i<=nener;i++)
	    xfrec[ib][i]=frec[ib][i]=0;

    h=1.0/(ymax-ymin);
    c1= 1.0 /(float) nener / h;
    c2= -0.5 /(float) nener / h +ymin;
					   /* segunda lectura para hacer */
					   /* un histograma              */

    if (argc==8)
    {
        Fnombres=fopen(nombres,"rt");
        if (Fnombres==NULL)
        {
            printf("\nEl fichero %s no existe\n\a",nombres);
            exit(0);
        }
    }


    Foutplt=fopen("analizado.plt","w");

    Fh=fopen("h.plt","w");

    titulo(nombre,nombre_case,&datos, n1, n2, l, nblo);
    fprintf(Foutplt,formato,nombre,nombre_case,
		    n_ener," -    -","MC it","","9P");

    numit=nblo*lblo*datos.itmax;
    itfrec=numit/500;
    if (itfrec==0)
	itfrec=1;
    ittot=0;
    for (n=n1;n<=n2;n++)
    {
      if (argc==8)
	fscanf(Fnombres,"%s\n",name);
      else
	sprintf(name,"%s%03d%s","OUT",n,".DAT");

      Finput=fopen(name,"rb");
      printf("%s\r",name);

      fread(&datos,sizeof(datos),1,Finput);


      for(idat=0;idat<jdat+1;idat++)
	fread(&v_dat[idat][0],4,(size_t) datos.itmax,Finput);

      fclose(Finput);
      ib=(n-n1)/lblo;
      
      for (it=0;it<datos.itmax;it++)
	{


	  y=v_dat[jdat][it];
	  x=(y-ymin)*h;
	  i=(int) (x*nener);
	  frec[ib][i]++;
	    
	  if (ittot%itfrec==0)
	    fprintf(Foutplt,"%ld %f\n",ittot,y);
	  ittot++;
	  xfrec[ib][i]+=v_dat[jdat][it];
  
	  
	  y=v_dat[jdat][it];
	  x=(y-ymin)*h;
	  i=(int) (x*nener);
	  frec[ib][i]++;
	  if (ittot%itfrec==0)
	    fprintf(Foutplt,"%ld %f\n",ittot,y);
	  ittot++;
	  xfrec[ib][i]+=v_dat[jdat][it];
        }
    
      
    }
    fprintf(Foutplt,"join linear\n");
    if (argc==8)
      fclose(Fnombres);
   
   
    /* dibuja histograma */
    titulo(nombre,nombre_case,&datos, n1, n2, l, nblo);
    fprintf(Foutplt,formato,nombre,nombre_case,
	    "n","",n_ener," -    -","9P");
    normalizacion=0;
    for(i=0;i<=nener;i++)
      for(j=0;j<nblo;j++)
	normalizacion+=(double) frec[j][i];
     normalizacion*=c1; 
    for (i=0;i<=nener;i++)
      {
	sumf=sumf2=0;
	for(ib=0;ib<nblo;ib++)
	  {
	    sumf+=frec[ib][i];
	    sumf2+=frec[ib][i]*frec[ib][i];
	  }
	/* normalizacion = 1;   sin normalizar */
	y=c1*i+c2;
	media=(double) sumf/ (double) nblo;
	disp=sqrt((sumf2/(double) nblo-media*media)/(double) (nblo-1))*
	  (double )nblo;
	hist[i]=sumf+disp;
	fprintf(Foutplt,"%10.6f %10.6f %10.6f\n",y,
		(double) sumf/normalizacion,disp/normalizacion);
	fprintf(Fh,"%10.6f %10.6f\n",y,
		(double) sumf/normalizacion);
      }
    fclose(Fh);

    fprintf(Foutplt,"bargraph shade 1\n");
    for(i=0;i<nener;i++)
      {
	y=c1*i+c2;
	fprintf(Foutplt,"%10.6f %10.6f\n",y,hist[i]/normalizacion);  
      }                                            /* normalizacion */
    fprintf(Foutplt,"join\n");
    
    delta= 2.0/(sigma*v);          /* maximo desplazamiento en FS */
    
    h=2.*delta/nbetas;
    y0=en;
    x0=coup;
    
    if (frac>0)
      /*	Comienza el ajuste de los picos a cubicas */
      {
	
	ndata=nener+1;
	for (i=0;i<ndata;i++)
	  {
	    e[i]=c1*i+c2;
	    for (ib=f[i]=0;ib<nblo;f[i]+=frec[ib++][i]);
	  }
	e0=e[ndata/2];
	vsol[0]=0;
	
	coupn=coup;
	
	coupn=datos_coupling[idat];

	for (j=0;j<ndata;j++) 
	  {
	    f[j]*=exp((e[j]-e0)*shift);
	    printf("\n f[%d]=%f",j,f[j]);
	  }
#ifdef MANUAL
	printf("introduce 4 numeros de 1 a 100 para el ajuste de los picos\n");
	scanf("%d %d %d %d",&used[0],&used[1],&used[2],&used[3]);
	for(i=0;i<4;i++)
	  used[i]=used[i]*ndata*0.01;	  
#endif
	for (i=0;i<3;i++)
	  {
	    
	    cubic_fit(e,f,frac,ndata,vsol,vb,used);
	    for(k=0;k<2;k++)
	      {
		xx=vsol[k];
		fr[k]=vb[4*k]+xx*(vb[4*k+1]+xx*(vb[4*k+2]+xx*vb[4*k+3]));
	      }
	    shift=log(fr[0]/fr[1])/(vsol[1]-vsol[0]);
	    coupn+=shift/v;
	    for (j=0;j<ndata;j++) f[j]*=exp((e[j]-e0)*shift);
	  }
      }
    else  /* jdat_c==3 */
      coupn=coup;
    cubic_fit(e,f,frac,ndata,vsol,vb,used);
    
    E1=vsol[0];
    E2=vsol[1];
    for(k=0;k<2;k++)
      {
	xx=vsol[k];
	fr[k]=vb[4*k]+xx*(vb[4*k+1]+xx*(vb[4*k+2]+xx*vb[4*k+3]));
      }
    FR1=fr[0];
    FR2=fr[1];
    suma=f[ndata-1]*(e[ndata-1]-e[ndata-2]);
    for (i=0;i<ndata-1;i++) suma += f[i]*(e[i+1]-e[i]);
    cte=1.0F/suma;
    
    fprintf(Foutplt,formato,nombre,nombre_case,
	    "f","",n_ener," -    -","4P");
    
    for (i=0;i<ndata;i++) fprintf(Foutplt,"%f %f\n",e[i],f[i]*cte);
    fprintf(Foutplt,"plot\n");
    
    for (k=0;k<2;k++)
      {
	for (i=used[2*k];i<used[2*k+1];i++)
	  {
	    xx=e[i];
	    yy=vb[4*k]+xx*(vb[4*k+1]+xx*(vb[4*k+2]+xx*vb[4*k+3]));
	    fprintf(Foutplt,"%f %f\n",xx,yy*cte);
	  }
	fprintf(Foutplt,"join\n");
      }
    
    sumL=sumLL=sumC=sumCC=0;
    for(k=0;k<2;k++)
      sumE[k]=sumEE[k]=0;
    for (ib=0;ib<nblo;ib++)
      {
	shift=(coupn-coup)*v;
	for (j=0;j<ndata;j++)
	  {
	    sum=0;
	    for (ibb=0;ibb<nblo;ibb++)
	      if (ibb!=ib) sum+=frec[ibb][j];
	    f[j]=sum*exp((e[j]-e0)*shift);
	  }
	cubic_fit(e,f,0.0,ndata,vsol,vb,used);
	for(k=0;k<2;k++)
	  {
	    xx=vsol[k];
	    fr[k]=vb[4*k]+xx*(vb[4*k+1]+xx*(vb[4*k+2]+xx*vb[4*k+3]));
	  }
	coup_b=log(fr[0]/fr[1])/(vsol[1]-vsol[0])/v;
	for (k=0;k<2;k++)
	  {
	    sumE[k]+=vsol[k];
	    sumEE[k]+=vsol[k]*vsol[k];
	  }
	xx=(vsol[1]-vsol[0]);
	sumL+=xx;
	sumLL+=xx*xx;
	sumC+=coup_b;
	sumCC+=coup_b*coup_b;
      }
    for (k=0;k<2;k++)
      {
	sumE[k]/= (double) nblo;
	E_err[k]=sqrt((sumEE[k]/(double) nblo-sumE[k]*sumE[k])* (double) (nblo-1));
      }
    sumL/=(double) nblo;
    sumLL/=(double) nblo;
    sumC/=(double) nblo;
    sumCC/=(double) nblo;
    coupn_err=sqrt((sumCC-sumC*sumC)*(double) (nblo-1));
    lat_err=sqrt((sumLL-sumL*sumL)*(double) (nblo-1));
    
    printf("\n E1=%6.4f+/-%6.4f; E2=%6.4f+/-%6.4f; Lat=%6.4f+/-%6.4f\n",
	   E1,E_err[0],E2,E_err[1],sumL,lat_err);
    printf("\n coupn =%8.6f +/- %8.6f\n",coupn,coupn_err);
    
    fprintf(Foutplt,"set order x dx y \n");
    fprintf(Foutplt,"set symbol 1P\n");
    fprintf(Foutplt,"%f %f %f\n",E1,E_err[0],FR1*cte);
    fprintf(Foutplt,"%f %f %f\n",E2,E_err[1],FR2*cte);
    fprintf(Foutplt,"plot \n");
    
    if (jdat_c!=3)
      fprintf(Foutplt,"title right font T 'Calor Latente=%6.4f+/-%6.4f'\n",
	      sumL,lat_err);
    else
      fprintf(Foutplt,"title right font T 'Magnetizacion=%6.4f+/-%6.4f'\n",
	      sumL,lat_err);
    


fclose(Foutplt);

}
    





