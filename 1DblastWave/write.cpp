#include "define.h"
#include "function.h"

void write(FILE *gp, int n, double *x, double *rho, double *rhou, double *E, int output_step)
{
	int i;
	double e[nn];
	double u[nn],p[nn];
	char fname[25];

	sprintf(fname,"%s%d%s","data/output_",n/output_step,".csv");
	gp = fopen(fname,"w");
	if(gp == NULL){
		printf("Can not open %s",fname);
		exit(1);
	}
	cal_up(rho,rhou,E,u,p);
	
	for(i=NGST;i!=NX+NGST;i++){
		e[i] = (E[i] - 0.5*rhou[i]*rhou[i]/rho[i])/rho[i];//e[i] = p[i]/(gamma-1.0)/rho[i];
		fprintf(gp,"%f,%f,%f,%f,%f,%f\n",x[i],rho[i],rhou[i],e[i],u[i],p[i]);
	}
	fclose(gp);
}
