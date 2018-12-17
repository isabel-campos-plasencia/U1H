#include "nrg.h"

extern short int u[];
extern int idel;
extern float del2;
extern float ebco[],mebco[];
extern int neigh_p[],neigh_m[];
extern int stap[];

extern int offset;

extern float good,tot;

void Staples(int j, int direccion)  /* calculo de las plaquetas incompletas */
{
    int ortog,offortog;
    int ja,jb,jc,jd;
    int nst;
    int vecino,offset;

    nst=0;
    vecino=neigh_p[direccion];
    jb=j+vecino;
    offset=direccion*V;
    for (ortog=0;ortog<4;ortog++)
	if (ortog!=direccion)
	{
	    offortog=ortog*V;
	    ja=j+neigh_p[ortog];                   
	    stap[nst++]= u[j+offortog]+u[ja+offset]-u[jb+offortog];
	    jc=j+neigh_m[ortog];
	    jd=jc+vecino;
	    stap[nst++]=-u[jc+offortog]+u[jc+offset]+u[jd+offortog];
	}
}

void Metropolis(int j, int direccion)      /*  algoritmo de Metropolis */
{
    int ihit,variacion;
    int uold,unew;
    int nst;
    float boltz_old,boltz_new;
    int offset;

    offset = direccion*V;
    uold=u[j+offset];
    boltz_old= ebco[Nestm1 & (stap[0]-uold)]
	      *ebco[Nestm1 & (stap[1]-uold)]
	      *ebco[Nestm1 & (stap[2]-uold)]
	      *ebco[Nestm1 & (stap[3]-uold)]
	      *ebco[Nestm1 & (stap[4]-uold)]
	      *ebco[Nestm1 & (stap[5]-uold)];

    for (ihit=0;ihit<nhit;ihit++)
    {
      tot++;
	variacion=( (int) (del2*FRANDOM) ) -idel;
	unew=uold+variacion;
	boltz_new=ebco[Nestm1 & (stap[0]-unew)]
		 *ebco[Nestm1 & (stap[1]-unew)]
		 *ebco[Nestm1 & (stap[2]-unew)]
		 *ebco[Nestm1 & (stap[3]-unew)]
		 *ebco[Nestm1 & (stap[4]-unew)]
		 *ebco[Nestm1 & (stap[5]-unew)];

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
}


