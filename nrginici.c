#include "nrg.h"

extern int x_p[],y_p[],z_p[],t_p[],
  x_m[],y_m[],z_m[],t_m[],v_p[],v_m[];

#ifdef HISTER
extern char dir[];
#endif
extern struct s_datos datos;


extern short int u[];

extern float coseno[],seno[],ebco[],mebco[];
extern float del2,delta,resultados[];
extern short int idel;
extern int normaener;
/*
extern int index0(int,int x[]);
extern void llenax_m(int, int, int);
extern int indexL(int, int x[]) ;
*/

/*int x[4];
 */
/*void fillcord(int, int, int, int cord[]);*/
int index0(int, int x[]);
void llenax_m(int , int, int );
int indexL(int, int x[]) ;
/*short int cordi[V][5];*/



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




int Direccionamientos(void)
{
  int x,y,z,t,cord[4],i,i1/*,index[V]*/,k,xin,xfin,yin,yfin,zin,zfin,tin,tfin;
  int l,i_0,cordl[4],i2,j,i_L,jx;
  FILE *fopen(),  *fp;







  for(i=0;i<V;i++)
    {                          
      x_m[i]=0;
      x_p[i]=0;
      y_m[i]=0; 
      y_p[i]=0;
      z_m[i]=0;
      z_p[i]=0;
      t_m[i]=0;
      t_p[i]=0;
      v_m[i]=0;
      v_p[i]=0;
    }
  /*  printf("%d",V);*/
  for(k=0;k<2;k++){
    /*    printf("\n k=%d ",k);*/
    for(i=0;i<D;i++){
      /*      printf("\n i=%d ",i);*/
      if(i==0&&k==0)xin=0;
      else xin=1;
      xfin=L;
      if(k==1&&i>0)xfin=Lm1;
      
      if(i<2&&k==0)yin=0;
      else yin=1;
      yfin=L;
      if(k==1&&i>1)yfin=Lm1;
      
      if(i<3&&k==0)zin=0;
      else zin=1;
      zfin=L;
      if(k==1&&i>2)zfin=Lm1;
      
      if(i<4&&k==0)tin=0;
      else tin=1;
      tfin=L;
      if(k==1&&i>3)tfin=Lm1;
	
      /*      printf(" 2");*/
      for (x=xin;x<xfin;x++){
	cord[0]=x;
	
	for(y=yin;y<yfin;y++){
	  cord[1]=y;
	  
	  for(z=zin;z<zfin;z++){
	    cord[2]=z;
	    
	    for(t=tin;t<tfin;t++){
	      cord[3]=t;
  /*   printf("xyzt = %1d%1d%1d%1d ",cord[0],cord[1],cord[2],cord[3]);*/
	      if(k){
		i1=indexL(i,cord);
		/*printf(" i1=%d ",i1);*/
/*fillcord(i1,i,k,cord);         */
		/*		index[i1]=i1;  */
		/*		i_0=0;         */
		i_L=-1;
		for(j=0;j<Dm1;j++){
		  if(j>=i&&cord[j]==Lm1&&i_L<0)i_L=j;
		  if(cord[j]>1){
		    cord[j]-=1;
		    i2=indexL(i,cord);
		    cord[j]+=1;
		    jx=(j<i)?j:j+1;
		    llenax_m(i1,i2,jx);
		  }
		  else{
		    if(j<i){
		      for(l=0;l<j;l++)
			cordl[l]=cord[l];
		      for(l=j;l<i;l++)
			cordl[l]=cord[l+1];
		      for(l=i;l<Dm1;l++)
			cordl[l]=cord[l];
		      cordl[i-1]=Lm1;
		      i2=index0(j,cordl);
		      llenax_m(i1,i2,j);
		      
		    }
		    else{
		      i_0=j+1;
		      for(l=0;l<i;l++)
			cordl[l]=cord[l];
		      for(l=i+1;l<i_0;l++)
			cordl[l]=cord[l-1];
		      for(l=i_0;l<Dm1;l++)
			cordl[l]=cord[l];
		      cordl[i]=Lm1;
		      i2=index0(i_0,cordl);
		      llenax_m(i1,i2,i_0);
		      
		    }
		  }
		  /*printf(" i2=%d ",i2);*/
		}
		/*printf(" 3");*/
		if(i_L>=0){
		  jx=i;
		  i_L+=1;
		  for(l=0;l<i;l++)
		    cordl[l]=cord[l];
		  for(l=i+1;l<i_L;l++)
		    cordl[l]=cord[l-1];
		  for(l=i_L;l<Dm1;l++)
		    cordl[l]=cord[l];
		  cordl[i]=Lm2;
		  i2=indexL(i_L,cordl);
		  llenax_m(i1,i2,jx);
		  
		}
		/*printf("i2=%d ",i2);*/
	      }
	      else{
		i1=index0(i,cord);
		/*printf("i1=%d ",i1);*/
/*		fillcord(i1,i,k,cord); */
		
		/*		index[i1]=i1;          */
		for(j=0;j<Dm1;j++){
		  if(cord[j]){
		    if(cord[j]==1&&j<i){
		      for(l=0;l<j;l++)
			cordl[l]=cord[l];
		      for(l=j;l<i-1;l++)
			cordl[l]=cord[l+1];
		      cordl[i-1]=0;
		      for(l=i;l<Dm1;l++)
			cordl[l]=cord[l];
		      i2=index0(j,cordl);
		      jx=(j<i)?j:j+1;
		      llenax_m(i1,i2,jx);
		    }
		    else{
		      cord[j]-=1;
		      i2=index0(i,cord);
		      cord[j]+=1;
		      jx=(j<i)?j:j+1;
		      llenax_m(i1,i2,jx);
		    }
		  }
		  /*printf("i2=%d ",i2);*/
		}
	      }
	    }

	    /*return(0);*/ /* no casca */
	  }
	  /*return(0);*/       
	}

      }

    }   /* for i */

  }
  /*printf(" 4");*/

return(0);  
}

int index0(int zero_i, int x[])
{
  int i,index2,index3,index4,Lm1_zero_i;
  

  for(i=0,index2=L_D;i<zero_i;i++){
    index2/=L;
    index2*=Lm1;
  }
  index3=0;
  Lm1_zero_i=1;
  if(zero_i){
    index3=x[zero_i-1]-1;
    Lm1_zero_i=Lm1;
    for(i=zero_i-1;i>0;i--){
      index3*=Lm1;
      index3+=x[i-1]-1;
      Lm1_zero_i*=Lm1;
    }
  }
  index4=0;
  if(D-zero_i>1){
    index4=x[D-2];
    for(i=D-2;i>zero_i;i--){
      index4*=L;
      index4+=x[i-1];
    }
  }
  return(L_D-index2+index3+index4*Lm1_zero_i);
}

int indexL(int L_i, int x[])
{
  int i,index1,index2,index3,index4,Lm2_L_i;
  
  /*Lm1=L-1;
    
    Lm2=Lm1-1;      */
  
  for(i=0,index2=Lm1_D;i<L_i;i++){
    index2/=Lm1;
    index2*=Lm2;
  }
  index3=0;
  Lm2_L_i=1;
  if(L_i){
    index3=x[L_i-1]-1;
    Lm2_L_i=Lm2;
    for(i=L_i-1;i>0;i--){
      index3*=Lm2;
      index3+=x[i-1]-1;
      Lm2_L_i*=Lm2;
    }
  }
  index4=0;
  if(D-L_i>1){
    index4=x[D-2]-1;
    for(i=D-2;i>L_i;i--){
      index4*=Lm1;
      index4+=x[i-1]-1;
    }
  }
  return(L_D-index2+index3+index4*Lm2_L_i);
}

void llenax_m(int i1, int i2, int jx)
{
  int di;
  
  /*if(i1>=V||i2>=V)printf("\n *********** Horror ******");*/
  di=i2-i1;
  switch (jx)  {
  case 0 :
    x_m[i1]=di;
    x_p[i2]=-di;
    break;
  case 1 :
    y_m[i1]=di;
    y_p[i2]=-di;
    break;
  case 2 :
    z_m[i1]=di;
    z_p[i2]=-di;
    break;
  case 3 :
    t_m[i1]=di;
    t_p[i2]=-di;
    break;
  case 4 :
    v_m[i1]=di;
    v_p[i2]=-di;
    break;
  }
}

void Inicializa(float beta, float gamma, int semilla, int flag)
{
  int i,x,y,z,t,dir,neigh_p[D],neigb_p[D],j,ja,jb;
  int igen,ortog,plaq,n4;
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
  idel=(short int) (delta/twopioverNest);
  del2= 2.0F * (float) idel + 1.0F;
  
  if (flag!=2)
    {
      igen=0;
      for(i=0;i<V;i++){
	
	if (flag == 0){
	  if(x_p[i])u[i]= (short int) (FRANDOM*rNest);
	  if(y_p[i])u[i+V]= (short int) (FRANDOM*rNest);
	  if(z_p[i])u[i+2*V]= (short int) (FRANDOM*rNest);
	  if(t_p[i])u[i+3*V]= (short int) (FRANDOM*rNest);
	  if(v_p[i])u[i+4*V]= (short int) (FRANDOM*rNest);
	}
	else{
	  if(x_p[i])u[i]= 0;
	  if(y_p[i])u[i+V]= 0;
	  if(z_p[i])u[i+2*V]= 0;
	  if(t_p[i])u[i+3*V]= 0;
	  if(v_p[i])u[i+4*V]= 0;
				}
      }
 
    }
  /*n4=0;*/
  normaener=0;
      for(j=0,plaq=0;j<V;j++){
	neigh_p[4]=v_p[j];
	neigh_p[3]=t_p[j];
	neigh_p[2]=z_p[j];
	neigh_p[1]=y_p[j];
	neigh_p[0]=x_p[j];
  for (dir=0;dir<D-1;dir++) /*   loop en direcciones de link    */
     {
	if(neigh_p[dir]){
	  jb=j+neigh_p[dir];
	  neigb_p[4]=v_p[jb];
	  neigb_p[3]=t_p[jb];
	  neigb_p[2]=z_p[jb];
	  neigb_p[1]=y_p[jb];
	  neigb_p[0]=x_p[jb];
	  for (ortog=dir+1;ortog<D;ortog++)
	    {
	      if(neigh_p[ortog]&&neigb_p[ortog]){
		normaener++;
		plaq++;}
	      
  	    }  /*for ortog*/
  	}   /*if neigh_p[dir]*/
	  /*plaq=0;*/
     } /*for j*/
      
      }  /*for dir*/


}

void Inicia_hister(float beta, float gamma, int semilla, int flag)
{
  int i,x,y,z,t,dir,neigh_p[D],neigb_p[D],j,ja,jb;
  int igen,ortog,plaq,n4;
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
  idel=(short int) (delta/twopioverNest);
  del2= 2.0F * (float) idel + 1.0F;
  

}


#ifdef HISTER
void escribe_hister(int i)
{

    FILE *Fhister;
    char name[100];
    int k;


    sprintf(name,"%s%s.dat",dir,"hister");

    Fhister=fopen(name,"a");    
    fprintf(Fhister,"%7.4f %7.4f ",datos.beta,datos.gamma);
    for (k=0;k<NOBS_HISTER;k++)
        fprintf(Fhister,"%f ",resultados[k]);
    fprintf(Fhister,"\n");    
    fclose(Fhister);   

    printf("(%5.3f,%5.3f),",datos.beta,datos.gamma);

}
#endif




