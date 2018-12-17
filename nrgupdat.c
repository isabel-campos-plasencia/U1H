#include "nrg.h"

/*extern int jrand,itp,itpf,jp,jpf,indmed;*/
/*  float en_i,en_f;*/
extern short int u[];
extern short int idel;
extern float del2;
extern float ebco[],mebco[];
extern int neigh_p[],neigh_m[];
extern short int stap[];
extern int x_p[],y_p[],z_p[],t_p[],
	   x_m[],y_m[],z_m[],t_m[],v_p[],v_m[];
extern  float coseno[];
extern void Medidas(void);
extern int offset;
extern float good;

int Staples(int j, int direccion)  /* calculo de las plaquetas incompletas */
{
  int ortog/*,offset*/,neigb_p[D],neigb_m[D];
  int ja,jb,jc,jd;
  int nst;
  int vecino,offortog;

  
  nst=0;
  vecino=neigh_p[direccion];
  jb=j+vecino;
  neigb_p[4]=v_p[jb];
  neigb_p[3]=t_p[jb];
  neigb_p[2]=z_p[jb];
  neigb_p[1]=y_p[jb];
  neigb_p[0]=x_p[jb];
  neigb_m[4]=v_m[jb];
  neigb_m[3]=t_m[jb];
  neigb_m[2]=z_m[jb];
  neigb_m[1]=y_m[jb];
  neigb_m[0]=x_m[jb];
  for (ortog=0;ortog<D;ortog++){
    if (ortog!=direccion)
      {
	offortog=ortog*V;
	if(neigh_p[ortog]&&neigb_p[ortog])
	  {
	    ja=j+neigh_p[ortog];  
	    /*    if((j+offortog)>=nlinks||(ja+offset)>=nlinks||(jb+offortog)>=nlinks)printf("\n ************* horrible *********** superior");*/
	    stap[nst++]= u[j+offortog]+u[ja+offset]-u[jb+offortog];
	    /*if(j==jrand&&itp==itpf&&jp==jpf)
{
printf("\n forstaps \n j=%d, dir=%d,offset = %d, -u[j,dir] =%d",j,direccion,offset,-u[j+offset]);
printf("\n ortog=%d,ja=%d,jb=%d,offortog=%d,u[j+ortog]=%d,u[ja]=%d,-u[jb]=%d",ortog,ja,jb, offortog,u[j+offortog],u[ja+offset],-u[jb+offortog]);
printf("\n nst=%d,stap=%d",nst-1,stap[nst-1]);
if(nst==1)
  {
Medidas();
getchar();
  }
}*/
	  }
	if(neigh_m[ortog]&&neigb_m[ortog])
	  {
	    jc=j+neigh_m[ortog];
	    jd=jb+neigb_m[ortog];
	    /*printf("\n offortog=%d",offortog);*/
	    /*  if((jc+offortog)>=nlinks||(jc+offset)>=nlinks||(jd+offortog)>=nlinks)
{
printf("\n ************* horrible *********** inferior");
printf("\n jc=%d offortog=%d offset=%d",jc,offortog,offset);
}*/
	    stap[nst++]=-u[jc+offortog]+u[jc+offset]+u[jd+offortog];
/*if(j==jrand&&itp==itpf&&jp==jpf)
{printf("\n backstaps \n jc=%d,jd=%d,-u[jc+ortog]=%d,u[jc]=%d,u[jd]=%d",jc,jd,-u[jc+offortog],u[jc+offset],u[jd+offortog]);
printf("\n nst=%d,stap=%d",nst-1,stap[nst-1]);
}*/
	  }
      }  /*if ortog*/
  }  /*for ortog*/
  return(nst);
}

void Metropolis(int j, int direccion, int nst)      /*  algoritmo de Metropolis */
{
  int ihit;
  short int variacion, uold,unew;
  int i;
  float boltz_old,boltz_new;


  
/*if(j==jrand&&itp==itpf&&jp==jpf)printf("\n nst en Metro = %d",nst);*/
  uold=u[j+offset];
  /*  en_i=0.;*/
  for(i=0,boltz_old=1;i<nst;i++)
    {
      boltz_old*= ebco[Nestm1 & (stap[i]-uold)];
      /*      en_i+=coseno[Nestm1 & (stap[i]-uold)];*/
    }
  /*indmed=0;*/
  /*Medidas();*/
  for (ihit=0;ihit<nhit;ihit++)
    {
      variacion=( (short int) (del2*FRANDOM) ) -idel;
      unew=uold+variacion;
      /*      en_f=0.0;*/
      for(i=0,boltz_new=1;i<nst;i++)
	{
	boltz_new*=ebco[Nestm1 & (stap[i]-unew)];
	/*  en_f+=coseno[Nestm1 & (stap[i]-unew)];*/
	}
      
      if (boltz_new < boltz_old)
	{
	  if (boltz_new > FRANDOM*boltz_old )
	    {
	      uold=unew;
	      boltz_old=boltz_new;
	      good++;
	    }
	}
      else
	{
	  uold=unew;
	  boltz_old=boltz_new;
	  good++;
	}
    }
  u[j+offset]=uold & Nestm1;
  /*indmed++;*/
  /*Medidas();*/
/*indmed++;*/
/*if(j==jrand&&itp==itpf&&jp==jpf)*/
  /*printf("\n ei =%f, ef =%f, ei-ef = %f",en_i,en_f,en_i-en_f);*/
}



