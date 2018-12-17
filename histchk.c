#include <stdlib.h>
#include <stdio.h>
#include <float.h>

void main(int argc, char *argv[])
{
    long int i,j,sum;
    float e;

    FILE *Finput;

    Finput=fopen(argv[1],"r");
    sum=0;
    while ((j=fscanf(Finput,"%f%ld",&e,&i))>0)
	sum+=i;
    printf(" total = %i \n",sum);
    fclose(Finput);
}