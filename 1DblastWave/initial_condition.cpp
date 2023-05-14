#include"define.h"
#include "function.h"

void initial(double *x, double dx, double *rho, double *rhou, double *E, double *rho_b, double *rhou_b, double *E_b)
{
	int i,ix;
	double x1,x2,u,p;
	
	for(i=0;i<nn;i++){
		if(i==0) x[i]  = xlow - (NGST * dx) + 0.5 * dx;
		else     x[i]  = x[i-1] + dx;
	}

//printf("x[NGST+NX/2]::  %f  %f %f  \n",x[NGST+NX/2-1],x[NGST+NX/2],x[NGST+NX/2+1]); 


if(ic==1) //sod
{
	x1 = 0.5;
	for(i=0;i<nn;i++){
		if(x[i]<=x1){
			rho[i]    = 1.0;
			u         = 0.0;
			p         = 1.0;
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];
		}else{
			rho[i]    = 0.125;
			u         = 0.0;
			p         = 0.1;
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];
		}
	}
}
else if(ic==2) //lax
{
	x1 = 0.5;
	for(i=0;i<nn;i++){
		if(x[i]<=x1){
			rho[i]    = 0.445;
			u         = 0.698;
			p         = 3.528;
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];
		}else{
			rho[i]    = 0.5;
			u         = 0.0;
			p         = 0.571;
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];
		}
	}

}
else if(ic==3)  //blast wave
{

		x1 = 0.1;
	x2 = 0.9;
	for(i=0;i<nn;i++){
		if(x[i]<x1){
			rho[i]    = 1.0;
			u         = 0.0;
			p         = 1000;
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];
		}
		else if(x[i]<x2)
		{
			rho[i]    = 1.0;
			u         = 0.0;
			p         = 0.01;
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];
		}
		else 
		{
			rho[i]    = 1.0;
			u         = 0.0;
			p         = 100.0;
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];

		}
	}
	
}
else if(ic==4) //density interaction
{

		x1 = -4.5;
	for(i=0;i<nn;i++){
		if(x[i]<=x1){
			rho[i]    = 1.515695;
			u         = 0.523346;
			p         = 1.805;
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];
		}else{
			rho[i]    = 1 + 0.1*sin(10*pi*x[i]);
			u         = 0.0;
			p         = 1.0;
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];
		}
	}



// end time 0.18
	// x1 = 0.1;
	// for(i=0;i<nn;i++)
	// {
	// 	if(x[i]<=x1)
	// 	{
	// 		rho[i]    = 3.857148;
	// 		u         = 2.629369;
	// 		p         = 10.333333;
	// 		rho_b[i]  = rho[i];
	// 		rhou[i]   = rho[i]*u;
	// 		rhou_b[i] = rhou[i];
	// 		E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
	// 		E_b[i]    = E[i];
	// 	}
	// 	else
	// 	{
	// 		rho[i]    = 1 + 0.2*sin(50*x[i]-25);
	// 		u         = 0.0;
	// 		p         = 1.0;
	// 		rho_b[i]  = rho[i];
	// 		rhou[i]   = rho[i]*u;
	// 		rhou_b[i] = rhou[i];
	// 		E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
	// 		E_b[i]    = E[i];
	// 	}
	// }
	
}
else if(ic==5) //Le Blanc 
{
		x1 = 3.0;
	for(i=0;i<nn;i++){
		if(x[i]<=x1){
			rho[i]    = 1.0;
			u         = 0.0;
			p         = 2.0/3.0*0.1;
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];
		}else{
			rho[i]    = 0.001;
			u         = 0.0;
			p         = 2.0/3.0*pow(10.0,-10);
			rho_b[i]  = rho[i];
			rhou[i]   = rho[i]*u;
			rhou_b[i] = rhou[i];
			E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
			E_b[i]    = E[i];
		}
	}
	
}
else if(ic==6) //Sedov
{

	for(i=0;i<nn;i++){

				rho[i]    = 1.0;
				u         = 0.0;
				p         = 1.0e-12;
				if(i==(NGST+NX/2))
				{
					//printf("x[i]::  %f  \n",x[i]); 
					p= (gamma-1.0)*rho[i]*3.2e6/dx;

				}
				rho_b[i]  = rho[i];
				rhou[i]   = rho[i]*u;
				rhou_b[i] = rhou[i];
				E[i]      = p / (gamma - 1.0) + 0.5 * rho[i] * u*u;
				E_b[i]    = E[i];
	}


	
}	


	

}
