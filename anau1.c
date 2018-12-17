
/*  PROGRAMA DE ANALISIS DEL MODELO SU(2) - HIGGS */


# include "anau1.h"

void muestra_datos(struct s_datos *dat)
{

	printf("\n El input de los archivos es: \n");
	printf("itmax  %d \n", dat->itmax);
	printf("mfresh  %d \n", dat->mfresh);
	printf("nbin  %d \n", dat->nbin);
	printf("itcut  %d \n", dat->itcut);
	printf("flag  %d \n", dat->flag);
	printf("seed  %d \n", dat->seed);
	printf("\n");
        printf("beta  %f \n", dat->beta);
        printf("gamma  %f \n", dat->gamma);
        printf("delta  %f \n", dat->delta);
	printf("\n");
}

void lee_argumentos( int argc, char *argv[],
		     struct s_datos *datos,
		     float *beta, float *gamma)
{
	int n,nbl,lbl;
	char nombres[100], name[100];
	FILE *Finput, *Fnombres;
	static struct s_datos datos_old;

	

	switch(argc)
	{
	 case 4: sscanf(argv[3],"%d", &nbl);
	 case 3: sscanf(argv[2],"%d", &n2);
		 sscanf(argv[1],"%d", &n1); break;
	 default: printf("ERROR: para lanzar introduce\n ");
		  printf(" ana%d (1er fich) (ult. fich)",L);
		  printf(" (n_bloques)\n ");
		  exit(0);
	 }
        nblo=nbl;

        iteraciones = (n2-n1 + 1);
        
	 if(n2-n1+1>maxfiles)
	 {
	    printf("MAXIMUM NUMBER OF FILES = %i \n", maxfiles);
	    exit(0);
	 }

	 printf("\n STARTING \n");
	 printf(" ANALYSIS OF OUT%03d.DAT - OUT%03d.DAT\n",n1,n2);
	 printf(" LATTICE SIZE %dx%dx%dx%d \n",L,L,L,L);
	 printf(" JACK-KNIFE BLOCKS %d\n ",nblo);

	 if((nblo)>maxbloque)
	 {
		printf("\n nblo  no puede superar %d \n",maxbloque);
		exit(0);
	 }

        
	 lbl=(n2-n1+1)/nbl;
	 if(lbl==0)
		exit(0);
	 lblo=lbl;
	 n1=n2-lblo*nblo+1;

	 for(n=n1;n<=n2;n++)
	 {
		sprintf(name,"%s%03d.DAT",nom_fich,n);
		if ((Finput=fopen(name,"rb"))==NULL)
		{
			printf("\n THIS FILE DOES NOT EXIST \n",name);
			exit(0);
		}

		fread(datos,sizeof(*datos),1,Finput);
	
		if(n==n1)
		{
		   memcpy(&datos_old,datos, sizeof(datos_old));
		}
		if(n>n1 && (datos->beta != datos_old.beta   ||
			    datos->gamma != datos_old.gamma   ||
			    datos->itmax != datos_old.itmax ||
			    datos->mfresh != datos_old.mfresh))
		{
			printf("NOT THE SAME INPUT FROM %d \n",n);
			exit (0);

		}
		fclose(Finput);

		if((ficheros[n-n1]= (char*)malloc(strlen(name)))==NULL)
		{
			printf("Fuera de memoria \n");
			exit(0);
		}
		strcpy(ficheros[n-n1],name);
	
          }

        printf("\n %d \n",datos->itmax);
        muestra_datos(datos);
        *beta = datos->beta;
        *gamma = datos->gamma;
        
}

void lee_datos(int n)
{
    int idat,idatnew,imag,j,i,k;
    FILE *Finput;
    struct s_datos dat;
     
        
    
    Finput=fopen(ficheros[n],"rb");
    fread(&dat,sizeof(dat),1,Finput);
         
    for(idat=0;idat<n_obs_medid;idat++)	     	     
        fread(&v_dat[idat][0],4,(size_t)dat.itmax,Finput);

    for(j=0;j<datos.itmax;j++)
      {
        v_dat[2][j] = v_dat[0][j]*v_dat[0][j];
        v_dat[3][j] = v_dat[2][j]*v_dat[2][j];

      }
        




                 
    fclose(Finput);
}

void Correlacion (int obs)
{
    int n,iter,itermax,it,i,id,maxdis;
    float ob[maxmed],cor[n_correl],corsum,sum1,sum2,med1[n_correl];
    float normal,med2[n_correl],r;
    FILE *Fcorrel,*FTau;
    char nombcorr[150],nomtau[100];
    double sum,mean,norm,tau;
    
    
    iter=0;
    for(n=n1;n<=n2;n++)
    {        
        lee_datos(n-n1);
        for(it=0;it<datos.itmax;it++)
            ob[iter++]=v_dat[obs][it];
    }
    itermax=iter;
    sprintf(nombcorr,"correl%d_%5.4f_%5.4f_%d.plt",obs,
            datos.beta,datos.gamma,L);
    Fcorrel=fopen(nombcorr,"w");
    fprintf(Fcorrel,"title left font T '% Corr.'\n");
    fprintf(Fcorrel,"title top font T 'Correlacion en tiempo de Montecarlo de %s'\n",cadena[obs]);
    fprintf(Fcorrel,"title bottom font T 'Iteraciones'\n");
    fprintf(Fcorrel, "set symbol 4N\n");
    
    sprintf(nomtau,"tau%d_%5.4f_%5.4f_%d.plt",obs,
            datos.beta,datos.gamma,L);
    FTau=fopen(nomtau,"w");
    fprintf(FTau,"title left font T '% tau'\n");
    fprintf(FTau,"title top font T 'Tau int, exp %s'\n",cadena[obs]);
    fprintf(FTau,"title bottom font T 'Iteraciones'\n");
    fprintf(FTau, "set symbol 4N\n");

    for(it=0;it<n_correl;it++)
    {
        corsum=0.;
        sum1=0.;
        sum2=0.;
        for(n=0;n<(itermax-it);n++)
        {
            corsum+=ob[n]*ob[n+it];
            sum1+=ob[n];
            sum2+=ob[n+it];
        }
        cor[it]=corsum/(itermax-it);
        med1[it]=sum1/(itermax-it);
        med2[it]=sum2/(itermax-it);
        
        cor[it]-=med1[it]*med2[it];
    }
    normal=cor[0];
    for(n=0;n<n_correl;n++)
    {
        cor[n]=cor[n]/normal;
        fprintf(Fcorrel,"%d\t%f\n",n,cor[n]);
    }
    
    fprintf(Fcorrel,"plot \n");
    fprintf(Fcorrel,"join \n");
    fclose(Fcorrel);

    /* CALCULO AUTOCONSISTENTE DE TAU */
    
    for(i=0;i<itermax;i++)
        sum +=ob[i];

    mean = sum / itermax;
    for(id=0;id<n_correl;id++)
    {
        sum =0;
        for(i=0;i<itermax;i++)
            sum += (ob[i] - mean)*(ob[i+id] -mean);
        
        if(!id) norm = sum;
        cor[id] = sum/norm;
        if( cor[id] <=0)
        {
            id++;
            break;
        }
    }
    maxdis = id;
    tau=0.5;
    for(id=1;id < maxdis; id++)
    {
        tau += cor[id];
        if(6*tau<id)
            break;
    }

    for(id=1;id < maxdis; id++)
    {
      r=cor[id-1]/cor[id];
      if(r>1)
	fprintf(FTau,"%d\t%lf\n",id,1.0/log(r));
    } 

    fprintf(FTau,"plot\n");
     for(id=1;id < maxdis; id++)
    {
      r=cor[id];
      if(r>0)
	fprintf(FTau,"%d\t%lf\n",id,(-1.0)*id/log(r));
    } 

    fprintf(FTau,"plot\n");
    fclose(FTau);


    printf("Tau integrado de %d = %f\n",obs,tau);
    
}  


void Histograma(int acoplo)
{
          int ib,j;
	  float frecu[n_inter+1],y;
	  FILE *Fhistog;
	  char nomhist[150];

	  sprintf(nomhist,"h%d_%5.4f_%5.4f_%d.plt",acoplo,
		  datos.beta,datos.gamma,L);
	  Fhistog=fopen(nomhist,"w");
	  fprintf(Fhistog,"title top fon T 'L=%d, b=%5.4f, g=%5.4,Distrib. de %d '\n",L,datos.beta,datos.gamma,acoplo);
	  fprintf(Fhistog,"title left font T 'frecuencia'\n");
	  fprintf(Fhistog,"title bottom font T 'Espectro de %d'",acoplo);
      
	  for(j=0;j<=n_inter;j++)
	    frecu[j]=0.;
	  for(j=0;j<=n_inter;j++)
	   for(ib=0;ib<nblo;ib++)
	    frecu[j]+=frec[ib][j];
	  for(j=0;j<=n_inter;j++)
	  {
	    y=c1*j+c2;
	    fprintf(Fhistog,"%f\t%f\n", y,frecu[j]);
	  }
	  fprintf(Fhistog,"bargraph solid fullwidth\n");
	  fclose(Fhistog);
}
Binder_Par(int nE)
{

  double binder,x,cumul;
  int j,imag,vol;
  FILE *FBINDER;
  char name1[100],name2[100],name3[100];
  double sm2,sm4,avm2_2,avm2_4,avm2_6,avm4_2;
  double sm,binder_error,tmp2,tmp4,tmp1,s1,der_error;


  vol=( (float) L*L*L*L);


  sprintf(name1,"E%dbener_%d.plt",nE,L);
  
  FBINDER=fopen(name1,"w");
  

  for(j=0;j<=nbetas;j++)
  {
          
     cumul = (double) (O_v[3][j] / (O_v[2][j]*O_v[2][j]));
      
      
      B[j] = (double) ( 1.0 - (double) (cumul/3.0) );
      
      
      /* calculo del error en el cumulante y en la derivada */
      
      /************************************************************/
      
      sm2 = O_err[2][j]*O_err[2][j];
      sm4 = O_err[3][j]*O_err[3][j];
      avm2_2 = O_v[2][j]*O_v[2][j];
      avm2_4 = avm2_2 * avm2_2;
      binder_error = ( (2*O_v[2][j]*sm2*O_v[3][j] + 
                        sm4*avm2_2) ) / (avm2_4) ;
      
      /************************************************************/  
      
      B_err[j] = sqrt(binder_error/3.0);    
      
      
      x = x0 - (delta*nbetas) + j*2*delta;
      
      fprintf(FBINDER,"%f\t%lf\t%lf\n",x,B[j],B_err[j]);  
      
      
  }
  fprintf(FBINDER, "join\n");
  fclose(FBINDER);
            
}


void FS(void)
{
	int ib,ibb,i,j,k,q,p,volu;
	double SumO[n_obs_FS],SumO2[n_obs_FS],SumderO[n_obs_FS],
               SumderO2[n_obs_FS],SumderlO[n_obs_FS],SumderlO2[n_obs_FS];
	double sumO[n_obs_FS],sumOe[n_obs_FS],sumx[n_obs_FS];
	double expo,expo_frec,O[n_obs_FS],derO[n_obs_FS],derlO[n_obs_FS];
	long int sumf;
  	double en, sum,sume,x,y,lO[n_obs_FS];
        char name[100];
  
    
  
        for(k=0;k<n_obs_FS;k++)
         for(ib=0;ib<=nblo;maxder[k][ib++]=0);

 	volu=L*L*L*L;
  

	for(j=0;j<=nbetas;j++)
	{
	 for(k=0;k<n_obs_FS;k++)
	 {
		SumO[k]=0.;
		SumO2[k]=0.;
		SumderO[k]=0.;
		SumderO2[k]=0.;
               
	 }

	 for(ib=0;ib<=nblo;ib++)
	 {
		sum=sume=0.;

		for(k=0;k<n_obs_FS;k++)
		 {
		  sumO[k]=0.;
		  sumOe[k]=0.;

		 }

		 x= x0 - (delta*nbetas) + j*2*delta;
               
                
		 for(i=0;i<n_inter;i++)
		 {
			sumf=0.;
			for(k=0;k<n_obs_FS;k++)
				sumx[k]=0.;
			y=c1*i+c2;
			for(ibb=0;ibb<nblo;ibb++)
				if(ib!=ibb)
				{
					sumf += frec[ibb][i];

					for(k=0;k<n_obs_FS;k++)
						sumx[k] += xfrec[k][ibb][i];
				}

			expo=exp((x-x0)*(y-ymed)*vol);
			expo_frec = expo*sumf;
			sum += expo_frec;
			sume += y*expo_frec;

			for(k=0;k<n_obs_FS;k++)
			{
				sumO[k] += sumx[k]*expo;
				sumOe[k] += sumx[k]*y*expo;
			}
		 }

		 en=sume/sum;
		 for(k=0;k<n_obs_FS;k++)
		 {
		     O[k]=sumO[k]/sum;
	             derO[k]=(sumOe[k]/sum - O[k]*en)*vol;
           	    
		 }

		
		 if(ib<nblo)
		  for(k=0;k<n_obs_FS;k++)
		  {
			SumO[k] += O[k];
			SumO2[k] += O[k]*O[k];
                        SumderO[k] += derO[k];
                        SumderO2[k] += derO[k]*derO[k];
                        
		   }
			
                 for(k=0;k<n_obs_FS;k++)
		    if(fabs(derO[k]) > fabs(maxder[k][ib]))
		    { 
                        
                       maxder[k][ib]=derO[k];
                       coup_maxder[k][ib]=x;
		    }

		 
	     } /* nblo */
		for(k=0;k<n_obs_FS;k++)
		{
			SumO[k]/=nblo;
		        SumderO[k]/=nblo;
			SumderlO[k]/=nblo;
                        O_v[k][j]=O[k];
		        derO_v[k][j]=derO[k];
			O_err[k][j]=
                                   sqrt(fabs(SumO2[k]/nblo-SumO[k]*SumO[k])*
                                   (double)(nblo-1));

                        derO_err[k][j]=
                                   sqrt(fabs(SumderO2[k]/nblo -
                                   SumderO[k]*SumderO[k])*(double)(nblo-1));

                        

		}


             
	} /* bucle en betas */

}


void Maxder(void)
{
	int ib,ibb,i,j,k,q,p,volu;
	double SumO[n_obs_FS],SumO2[n_obs_FS],SumderO[n_obs_FS];
        double SumderO2[n_obs_FS];
	double sumO[n_obs_FS],sumOe[n_obs_FS],sumx[n_obs_FS];
	double expo,expo_frec,O[n_obs_FS],derO[n_obs_FS];
        long double Sum, Sum2,Sumder,Sumder2;
	long int sumf;
	double en, sum,sume,x,y;


	volu=L*L*L*L/2;
	 for(ib=0;ib<=nblo;ib++)
	 {

	  for(k=0;k<n_obs_FS;k++)
	  {

		  x=coup_maxder[k][ib]-h/2-h/20;
		  for(j=-10;j<=10;j++)
		  {
			x+=h/20; 

			if(k<n_obs_FS)
			{

			  sumO[k]=0.;
			  sumOe[k]=0.;
			  sum=sume=0.;


			 for(i=0;i<n_inter;i++)
			 {
				sumf=0.;
			      
				sumx[k]=0.;

				for(ibb=0;ibb<nblo;ibb++)
					if(ib!=ibb)
					{
					sumf += frec[ibb][i];
					sumx[k] += xfrec[k][ibb][i];
					}

				expo=exp((x-x0)*(y-ymed)*vol);
				expo_frec = expo*sumf;
				sum += expo_frec;
				sume += y*expo_frec;

	
	
					sumO[k] += sumx[k]*expo;
					sumOe[k] += sumx[k]*y*expo;
				
			 }

			 en=sume/sum;
		      
			 
			 O[k]=sumO[k]/sum;

/* esto es la derivada */

		         derO[k]=(sumOe[k]/sum-O[k]*en)*vol;
			 
			}
			
		   
			if((derO[k])>(maxder[k][ib]))
			{
				maxder[k][ib]=derO[k];
				coup_maxder[k][ib]=x;
			}
	  }/* bucle en betas */
  

	 } /* bucle en observable  */
    
	} /* bucle en bloques */

	for(k=0;k<n_obs_FS;k++)
	{
		Sum=Sum2=Sumder=Sumder2=0.;
		for(ib=0;ib<nblo;ib++)
		{
			Sum += coup_maxder[k][ib];
			Sum2 += coup_maxder[k][ib]*coup_maxder[k][ib];
			Sumder += maxder[k][ib];
			Sumder2 += maxder[k][ib]*maxder[k][ib];
		}
		Sum/=nblo;
		err_coup[k]=sqrt(fabs(Sum2/nblo-Sum*Sum)*(double)(nblo-1));
		Sumder/=nblo;
		err_maxder[k]=sqrt(fabs(Sumder2/nblo-Sumder*Sumder)*(double)(nblo-1));
	}
}




void titulo (void)
{
  int i;
  sprintf(nombre,"title top font T 'L=%d, b=%f, g=%f, #MC= %dx%d Nbl=%d, Lblo=%d, %s'",L,datos.beta,datos.gamma,datos.itmax,datos.mfresh,nblo, lblo, nom);


}

void Dibuja_resultados (int v)
{
	int j;
	double x;
	char formato[]= "new plot \n"
	      		"%s \n"
		        "title left font T '%s' \n"
			"title right font T '%s' \n"
			"case               '%s' \n"
			"set order x y dy \n"
			"set labels left font T \n"
		        "set labels bottom font T \n"
			"set symbol %s \n";


	/* OBSERVABLE  */

	fprintf(Foutplt,formato,nombre,"O",n_coup,"G- -","9N");

	for(j=0;j<=nbetas;j++)
	{
		x=x0- (nbetas*delta) + (2*j*delta);
		fprintf(Foutplt,"%10f  %10f\n",x,O_v[v][j]);
	}
	fprintf(Foutplt,"join \n");

	for(j=0;j<=nbetas;j++)
	{

		x=x0-(nbetas*delta) + (2*j*delta);
		fprintf(Foutplt,"%10f  %10f  %10f\n",x,O_v[v][j], O_err[v][j]);
	}
	fprintf(Foutplt,"plot\n");
	fprintf(Foutplt,"set symbol 5N \n");
	fprintf(Foutplt, " %10f %10f %10f\n", x0,O_v[v][nbetas/2],
		O_err[v][nbetas/2]);
	fprintf(Foutplt, "plot \n");

	/* DERIVADA  */
	fprintf(Foutplt, formato, nombre,n_deri,n_coup,"G- -","9N");
    
	for(j=0;j<=nbetas;j++)
	{
	  x=x0- (nbetas*delta) + (2*j*delta);
	  fprintf(Foutplt,"%10f  %10f\n",x,derO_v[v][j]);
	}
	fprintf(Foutplt,"join \n");

	for(j=0;j<=nbetas;j++)
	{

	    x=x0-(nbetas*delta) + (2*j*delta);
	    fprintf(Foutplt,"%10f  %10f  %10f\n",
		    x,derO_v[v][j], derO_err[v][j]);
	}
	
       
        fprintf(Foutplt,"plot\n");
	fprintf(Foutplt,"set symbol 3N \n");
	fprintf(Foutplt,"set order x dx y dy \n");
	fprintf(Foutplt,"%10f %10f %10f %10f\n",coup_maxder[v][nblo],
		err_coup[v],maxder[v][nblo],err_maxder[v]);
	fprintf(Foutplt,"set order x y dy \n");
	fprintf(Foutplt, " %10f %10f %10f \n", x0,derO_v[v][nbetas/2],
		derO_err[v][nbetas/2]);
	fprintf(Foutplt, "plot \n");

	fprintf(Foutplt,"title bottom font T' C=%8.6f+/-%8.6f; D=%8f+/-%8f'\n",
		coup_maxder[v][nblo],err_coup[v],
		maxder[v][nblo],err_maxder[v]);


}

void Evolucion(void)
{       
          FILE *Fevol[n_obs_plot];
          int n,it,obs;
          char nombreevol[n_obs_plot][150];

          for(obs=0;obs<n_obs_medid;obs++)   /* E y parametros de orden */
          {  
            sprintf(nombreevol[obs],"evol%d.plt",obs);
            if((Fevol[obs]=fopen(nombreevol[obs],"w"))==NULL)
              {
                puts("No puedo abrir el archivo de resultados\n");
                exit(1);
              }
          }

	  medidas = 5;

          for(n=n1;n<=n2;n++)
          {
           lee_datos(n-n1);

           for(it=0;it<datos.itmax;it++)
           {      
              if((it%medidas) == 0)
              {
               for(obs=0;obs<n_obs_medid;obs++)
                {
                  fprintf(Fevol[obs],"%d\t%lf\n",
			  ((n-n1)*datos.itmax+it),v_dat[obs][it]);
                }
              }
            }
          }

          for(obs=0;obs< n_obs_medid;obs++)
            {
              fprintf(Fevol[obs],"join\n");
              fclose(Fevol[obs]);
            }
}         





void main(int argc,char *argv[])
{

    char nombreresult[150],n_der[3][10],n_cou[3][6];
    
    char nomeigen[50],nombrefichero[150];

    int  n_vo[3],j,n,it,ib,t,i,p,Num_Ener,maxi,maxj;
    float datos_beta[3],sum,sumen,sumen2,ymin,ymax,beta,gamma,coup;
    float en,cv,sigma,y,aut[4],max,s_prom[4],promedio[4],desviacion[4];
    FILE *Foutput;
    int k,mm,kk,itt;
    
    
    strcpy(n_cou[0],"beta");
    strcpy(n_cou[1],"gamma");
    
    
    lee_argumentos(argc,argv,&datos,&beta,&gamma);
    
    n_vo[0]=6.0;
    n_vo[1]=6.0;            /* CHEQUEAR !! */


    datos_beta[0]=datos.beta;
    datos_beta[1]=datos.gamma;

    
    printf("Los observables analizados son:\n");
    for(j=0;j< (n_obs_medid);j++)
        printf("          %d\t%s\n",j,cadena[j]);
    
    Correlacion(0);
    Correlacion(1);

    
    Evolucion();
   
    for(d=0;d<2;d++) 
    {

        sprintf(nombreresult,"%dresult_%5.4f_%5.4f_%d.plt",
                d,beta,gamma,L);
        if((Foutplt=fopen(nombreresult,"w"))==NULL)
        {
            puts("No puedo abrir el archivo de resultados\n");
            exit(1);
        }

	      
        strcpy(n_coup,n_cou[d]);
        n_vol=n_vo[d];
        coup=datos_beta[d];
        vol=(long int)n_vol*L*L*L*L;  
        
        
        sum=sumen=sumen2=0.;
        ymin=10.;               /* Valores maximo y minimo de E */
        ymax=-10.0;

        for(n=n1;n<=n2;n++)
        {
            lee_datos(n-n1);
            for(it=0;it<datos.itmax;it++)
            {
                y=v_dat[d][it];
                sum++;
                sumen += y;
                sumen2 += y*y;
                ymin=(y<ymin)?y:ymin;
                ymax=(y>ymax)?y:ymax;				
            }
        }
        en=sumen/sum;                  /* Dispersion de la gaussiana */
        cv=fabs(sumen2/sum-en*en);     /* para hacer FS              */
        sigma=sqrt(cv);
        
        for(n=0;n<n_obs_FS;n++)
            for(ib=0;ib<nblo;ib++)
                for(j=0;j<=n_inter;j++)
                    xfrec[n][ib][j]=frec[ib][j]=0.;
        
        h=1.0/(ymax-ymin);
        c1=1.0/(float)n_inter/h;
        c2=-0.5/(float)n_inter/h+ymin;
        
        for(n=n1;n<=n2;n++)
            
        {                            /* Clasificacion en intervalos */
            ib=(n-n1)/lblo;      /* de los Observ. */
            lee_datos(n-n1);
            for(it=0;it<datos.itmax;it++)
            {

                y=v_dat[d][it];
                i=(int)((y-ymin)*h*n_inter);
                
                frec[ib][i]++;
                for(t=0;t<n_obs_FS;t++)
                    xfrec[t][ib][i]+=v_dat[t][it];
		
            }     
        }
        delta=0.00001;
        
        ymed=en;
        x0=coup;
        Histograma(d);
        FS();
	Maxder();  

 
      
	if(d==0) Binder_Par(0);   
      
      
        

        for(n=0;n<(n_obs_FS );n++)
	{
            
            sprintf(nom," O= %s ",cadena[n]);
            titulo();
            Dibuja_resultados(n);
	}
 

        fclose(Foutplt);
    }

    
    
    for(n=0;n< (n_obs_FS  );n++)
      printf("\n %d. %s= %8.6f +/- %8.6f \n",n,cadena[n],
	     O_v[n][nbetas/2],O_err[n][nbetas/2]);
    
   



}          /*  END */
	












