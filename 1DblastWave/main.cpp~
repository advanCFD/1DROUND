#include "define.h"
#include "function.h"

int main(void)
{
	double t;
	int i,n;
	double rho[nn];  // ���x
	double rhou[nn]; // �^����
	double E[nn];    // �G�l���M�[
	double rho_b[nn],rhou_b[nn],E_b[nn]; // 1step�O
	double x[nn];    // �Z����x���W
	double dt,dx;
	double CFL;
	double total_t;
	int output_step;
	double xe[5000],fe[5000],pe[5000],ue[5000];
//	FILE *gp, *fp, *ep;
	
	CFL = 0.2;
	n = 0;
	t = 0.0;
	dt = 0.0;
	dx = (xup-xlow)/NX;
	output_step = 0;
	total_t = 0.038;
	output_step = 10;
	initial(x,dx,rho,rhou,E,rho_b,rhou_b,E_b);
//	fp = popen("gnuplot -persist", "w");

	FILE *write1;
	write1=fopen("1D_Euler.dat","w");
	fclose(write1);

	//Output(x,rho,rhou,E,0);
	//exit(0.0);
	
	while(t < total_t){ printf("t: %f  \n",t);
		for(i=0;i<nn;i++){
			rho_b[i] = rho[i];
			rhou_b[i] = rhou[i];
			E_b[i] = E[i];
		}
		if(time_in == 1) EF(n,&dt,dx,CFL,rho,rhou,E,rho_b,rhou_b,E_b,t,total_t);
		else if(time_in == 2) SSP_RK3(n,&dt,dx,CFL,rho,rhou,E,rho_b,rhou_b,E_b,t,total_t);
		else if(time_in == 3) RK5(n,&dt,dx,CFL,rho,rhou,E,rho_b,rhou_b,E_b,t,total_t);
		t = t + dt;
		if(n%output_step == 0) printf("time step : %f / %f dt = %f \n",t,total_t,dt);
		//if(n%output_step == 0) write(gp,n,x,rho,rhou,E,output_step);
		//graph(fp,x,rho);
		n++;
		//Output(x,rho,rhou,E,n);
	}
	printf("time step : %f / %f dt = %f \n",t,total_t,dt);
	/*
	graph(fp,x,rho);
	write(gp,n,x,rho,rhou,E,output_step);
	fprintf(fp, "exit\n");                        // gnuplot�̏I��
	pclose(fp);
	*/
	Output(x,rho,rhou,E,1);
}


void Output(double *x, double *rho, double *rhou, double *E, int step){
	FILE *write;
	
	write=fopen("1D_Euler.dat","a");
	fprintf(write,"%s","Title=\"ROUND RESULT\"\n");
	fprintf(write,"%s","VARIBLES=\"X\",\"Density\",\"Velocity\",\"Pressure\",\"Internal Energy\"\n");
	fprintf(write,"%s","ZONE T=\"STD: Geom\"\n");
	fprintf(write,"%s %d %s","STRANDID=",step+1,"\n");	
	fprintf(write,"%s %d %s","SOLUTIONTIME=",step,"\n");	
	fprintf(write,"%s %d %s","I=",NX,"\n");	
	fprintf(write,"%s ","F=POINT\n");	

	double u[nn],p[nn];
	double e[nn];
	cal_up(rho,rhou,E,u,p);
	for (int i=NGST;i!=NX+NGST;++i){

		e[i] = (E[i] - 0.5*rhou[i]*rhou[i]/rho[i])/rho[i];
		//fprintf(write,"%f %f %f %f %f %f%s",x[i],rho[i],rhou[i],e[i],u[i],p[i],"\n");
		fprintf(write,"%f %f %f %f %f %f%s",x[i],rho[i],u[i],p[i],rhou[i],e[i],"\n");
	}

	fclose(write);

}

