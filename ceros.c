
/* Dada la distribucion de probabilidad calcula los ceros
   de la funcion de particion ESFERA  */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

#define nbetas 100
#define ninter_eta 50
#define ninter_omega 200


float PE[nbetas],E[nbetas];

float delta_e,delta_w,beta,gamma;


int N;  /* lado del hipercubo 5D */


void lee_input(void)
{

    FILE *finput;

    finput=fopen("datos","r");
    if(finput==NULL)
    {
        printf("No existe datos \n");
        exit(1);
    }    
    fscanf(finput,"%f",&delta_e);
    fscanf(finput,"%f",&delta_w);
    fscanf(finput,"%d",&N);
    fscanf(finput,"%f",&beta);           /* beta de la simulacion */
    fscanf(finput,"%f",&gamma);
    printf("Calculo de la funcion en el intervalo: \n");
    printf("beta = %f\t %f\n",  beta-delta_e,beta+delta_e);
    printf("Imag  = %f\n", delta_w);
    printf("con acoplamientos: \n");
    printf("beta de la simulacion = %f\n",beta);
    printf("gamma de la simulacion = %f\n",gamma);
    printf("N = %d\n",N);

  }


void lee_distribucion(void)
{

    FILE *fdistrib;

    int i,j;

    fdistrib=fopen("h.plt","r");
    if(fdistrib==NULL)
    {
        printf("No esta la distribucion de Probabilidad \n");
        exit(1);
    }    

    for(i=0;i<nbetas;i++)
	fscanf(fdistrib,"%f\t%f",&E[i],&PE[i]);
    
}




void main()
{      
   
    FILE *fZ;
    
    int i,j,k,V;

    float eta,omega,pasose,pasosw,eta_minimo,omega_minimo;
    double c,s,denom,funcion,minimo;
   
    float EM;

    lee_input();

    /*  LA DIFERENCIA ESTA AL CALCULAR EL NUMERO DE PLAQUETAS  */

    V = 60*(N-1)*(N-1)*(N-1)*(N-1) - 20*(N-1)*(N-1);

    lee_distribucion();

    EM=0.0;

    for(j=0;j<nbetas;j++)
      EM += E[j];

    EM /= nbetas;


    pasose = (delta_e)/ninter_eta;
    pasosw = (delta_w)/ninter_omega;


    fZ=fopen("Z.plt","w");    
    fprintf(fZ,"set order dummy x y \n");
    fprintf(fZ,"set limit y -0.05 1 \n");

    eta_minimo = omega_minimo = 0.0;
    minimo = 1;

    /* Se calcula la prediccion de FS en beta'= eta + i*omega */

    for(i=1;i<2*ninter_eta;i++)
      {
	eta = (beta - ninter_eta*pasose) + i*pasose;

	denom = 0.0;

        for(j=0;j<nbetas;j++)
	  denom += PE[j]*exp(-(eta - beta)*V*(E[j]-EM));

	denom *= denom;

	for(k=1;k<ninter_omega;k++)
	  {
	    omega = k*pasosw;

	    c = s = 0.0;

	    for(j=0;j<nbetas;j++)
	      {
		c += PE[j]*exp(-(eta - beta)*V*(E[j]-EM))*cos(V*E[j]*omega);
		s += PE[j]*exp(-(eta - beta)*V*(E[j]-EM))*sin(V*E[j]*omega);
	    
	      }
	    c *= c;
	    s *= s;
	    
	    funcion = (c + s)/denom;

	    if(minimo > funcion) 
	      {
		minimo = funcion;
		eta_minimo = eta;
		omega_minimo = omega;
	      }

	    fprintf(fZ,"%f\t%f\t%f\n",eta,omega,funcion);

	  }

	fprintf(fZ,"plot\n"); 
       
      }


    printf("\n\n parte Imaginaria del minimo = %f\n",omega_minimo);
    printf("\n\n parte Real del minimo = %f\n",eta_minimo);
    printf("\n\n Valor de la funcion en el minimo = %f\n",minimo);
    fclose(fZ);


  }


