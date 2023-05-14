#include "define.h"
#include "function.h"

void period(double *f)
{
	int i;
	for(i=0;i!=NGST;++i){
		f[i]=f[i+NX];
		f[i+NX+NGST]=f[i+NGST];
	}
}


void wall(double *rho, double *rhou, double *E)
{
	int i;
	for(i=0;i<NGST;i++){
		rho[i]  = rho[NGST];
		rhou[i] = -1*rhou[NGST];//-1*rhou[NGST];
		E[i] = E[NGST];
	}
	
	for(i=NX+NGST;i<nn;i++){
		rho[i] = rho[NX+NGST-1];
		rhou[i] = -1*rhou[NX+NGST-1];//-1*rhou[NX+NGST-1];
		E[i] = E[NX+NGST-1];
		//printf("%3d rho = %f rhou = %f E = %f\n",i,rho[i],rhou[i],E[i]);
	}
}


void neumann(double *f)
{
	int i;
	for(i=0;i<NGST;i++) f[i] = f[NGST];
	for(i=NX+NGST;i<nn;i++) f[i] = f[NX+NGST-1];
}
