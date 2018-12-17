#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
/*#define DEBUG*/
/* #define HISTER */

#define L         16
#define Lm1			(L-1)
#define Lm2			(L-2)
#define maxglob    L
#define D			5
#define	Dm1		(D-1)
#define L_D 		(L*L*L*L*L)
#define Lm1_D		(Lm1*Lm1*Lm1*Lm1*Lm1)
#define Nest       1024        /* U(1) se aproxima por Z(Nest) */

#define nhit       3          /* Numero de hits de los c. de gauge         */

#define maxit      200        /* Numero de iteraciones maximo por bin      */

#define NOBS_HISTER 2			

#define  V         (L_D-Lm2*Lm2*Lm2*Lm2*Lm2)  /* numero de sites en cada transputer        */
/*#define  Voff		V-L^(D-1)-(L-2)^(D-1)  no usar               */
#define  nlinks    (D*V)


#define L3 (L*L*L)
#define L2 (L*L)

#define twopi 6.283185307

#define NormRAN (1.0F/( (float) RAND_MAX+1.0F))

#define  RAN() ( (float) rand() * NormRAN )


#define FNORM   (2.3283063671E-10F)
#define RANDOM  ( (ira[ip++]=ira[ip1++]+ira[ip2++]) ^ira[ip3++] )
#define FRANDOM (FNORM*RANDOM)


#ifdef MAIN
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
#else
extern unsigned char ip,ip1,ip2,ip3;
extern unsigned ira[];
#endif

#define randmax NormRAN



#define  Nestm1  (Nest-1)

#define twopioverNest ( twopi/((double) Nest))
#define Nestovertwopi ( ((float) Nest)/twopi)
/*#define Normaener  ( (float) (1.0/(double) (V*6)))         */
#define Normapoly  ( (float) (1.0/(double) (4*L3)))

#define rNest      ( (float) Nest)

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000L
#endif

#define n_obs 4

struct s_datos
{
   int itmax,           /* Numero de medidas por bloque                  */
       mesfr,           /* frecuencia de las medidas                     */
       nbin,            /* numero de bloque                              */
       itcut,           /* proximo bloque a calcular                     */
       flag,            /* conf de partida: 0(random), 1(fria),2(backup) */
       seed;            /* semilla random                                */
    float beta,         /* acoplamiento de los campos de gauge           */
          gamma,        /* acoplamiento plaquetas al cuadrado            */
	  delta        /* salto de Metropolis para los campos de gauge  */
#ifdef HISTER
        ,dbeta,        /* variacion de beta para la histeresis   */
         dgamma          /* variacion de gamma para la histeresis      */
#endif
             ;   
};

#define NDATOS_INT   6     /* numero de campos int   en s_datos */
#ifdef HISTER
#define NDATOS_FLOAT 5     /* numero de campos float en s_datos */
#else
#define NDATOS_FLOAT 3     /* numero de campos float en s_datos */
#endif









