
/* PROGRAMA DE ANALISIS PARA U(1) CON DOS COUPLINGS EN LA ESFERA */
		     /*    HEADER	*/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>


typedef struct{float a0,a1,a2,a3;} o4v;

 

# define L              18
# define Lm1            (L-1)
# define Lm2            (L-2)
# define V              (L*L*L*L*L-Lm2*Lm2*Lm2*Lm2*Lm2)
# define N_Pl (Lm1*(4+Lm2*(19+Lm2*(35+Lm2*(20-Lm1)-2*Lm1*Lm1))+Lm1*(4+L*(3+L*(2+L))+2*Lm1*(6+L*(3+L)+9*Lm1))))
/*# define N_Pl           5040*/ /* Valido para L=4 */
/*#define N_Pl            16800*/ /* Nueva Medida para L=4*/
# define maxit 		5000
# define maxbloque 	 900
# define nom_fich	 "OUT"
# define n_obs_medid     2 /* estos se miden directamente del output */    
# define n_obs_FS	 4 /* para el cumulante */
# define maxfiles	 5000
# define nbetas  	 80
# define n_inter	 100
# define n_correl	 100
# define maxmed		 200000
# define n_obs_plot      2
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
        a[k][l]=h+s*(g-h*tau);



float v_dat[n_obs_FS][maxit],raw_dat[n_obs_medid][maxit];
long int frec[maxbloque][n_inter+1];
long int vol=N_Pl;
double xfrec[n_obs_FS][maxbloque][n_inter+1],SumO[n_obs_FS];
double coup_maxder[n_obs_FS][maxbloque+1], err_coup[n_obs_FS];
double maxder[n_obs_FS][maxbloque+1], err_maxder[n_obs_FS];
double O_v[n_obs_FS][nbetas+1], O_err[n_obs_FS][nbetas+1];
double derO_v[n_obs_FS][nbetas+1], derO_err[n_obs_FS][nbetas+1];
double ymed,x0,h,c1,c2,delta,derlO_v[n_obs_FS][nbetas+1],
        derlO_err[n_obs_FS][nbetas+1];
char   *ficheros[maxfiles],n_coup[10],n_deri[10], nombre[150],nom[100];
int    n_vol,permutar,d,medidas,iteraciones,m;
double B[nbetas+1],B_err[nbetas+1];


FILE *Foutplt;

int nblo,lblo,n1,n2;


struct s_datos
{
 int	itmax,
		mfresh,
		  nbin,
		itcut,
		flag,
		seed;
 float  beta,
   gamma,
 delta;
 };

 struct s_datos datos;
 char cadena[n_obs_FS][100]={    
     "Energia plaquetas",
     "Energia plaquetas cuadrado"};








