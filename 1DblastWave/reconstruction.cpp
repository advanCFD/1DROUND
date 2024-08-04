#include "define.h"
#include "function.h"


double WENO(double w1, double w2, double w3, double w4, double w5)
{
	double omega0,omega1,omega2,omega;
	double S0=w1/3.0-7*w2/6.0+11.0*w3/6.0;
	double S1=-w2/6.0+5*w3/6.0+w4/3.0;
	double S2=w3/3.0+5*w4/6.0-w5/6.0;
	
	double beta0=13.0*pow(w1-2*w2+w3,2.0)/12.0+pow(w1-4*w2+3*w3,2.0)/4.0;
	double beta1=13.0*pow(w2-2*w3+w4,2.0)/12.0+pow(w2-w4,2.0)/4.0;
	double beta2=13.0*pow(w3-2*w4+w5,2.0)/12.0+pow(3*w3-4*w4+w5,2.0)/4.0;
	
	if(WENO_Method==0){
		omega0=0.1/pow(epsilon+beta0,2.0);
		omega1=0.6/pow(epsilon+beta1,2.0);
		omega2=0.3/pow(epsilon+beta2,2.0);
	}
	if(WENO_Method==1){
		double Tao5=fabs(beta2-beta0);
		omega0=0.1*(1.0+(Tao5/(beta0+epsilon1)));
		omega1=0.6*(1.0+(Tao5/(beta1+epsilon1)));
		omega2=0.3*(1.0+(Tao5/(beta2+epsilon1)));
	}
	omega=omega0+omega1+omega2;
	omega0=omega0/omega;
	omega1=omega1/omega;
	omega2=omega2/omega;
	
	return omega0*S0+omega1*S1+omega2*S2;
}

double ROUND(double u1, double u2, double u3)
{

  double zeps  = 1.0e-20;

  double z0=(u2-u1+zeps)/(u3-u1+zeps);
  
  double g=z0;

  if(ROUND_Method==0){
  
    if(z0>0.0 && z0<=0.5){

    double a=1100.0;
    double w=1.0/(1.0+a*square(square(z0-0.5)))/(1.0+a*square(square(z0-0.5)));
    double gp=(1.0/3.0+5.0/6.0*z0)*w+2.0*z0*(1.0-w);
    g=min(gp,2.0*z0);

    }

    if(z0>0.5 && z0<1.0){

    double a=800.0;
    double w=1.0/(1.0+a*square(square(z0-0.5)))/(1.0+a*square(square(z0-0.5)));
    double gp=(1.0/3.0+5.0/6.0*z0)*w+(0.15*z0+0.85)*(1.0-w);
    g=min(gp,0.15*z0+0.85);

    }
  
  }
  
  
  if(ROUND_Method==1){
  
  double w1=1.0/square(1.0+25.0*z0*z0);
  
  double w2=1.0/square(1.0+25.0*(1.0-z0)*(1.0-z0));
  
  g=(5.0/6.0+2.0/3.0*w1-1.0/3.0*w2)*z0+1.0/3.0-1.0/3.0*w1+1.0/6.0*w2;
  
  
  }  

  if(ROUND_Method==2){
  
    if(z0<=0.0){
    double a=-0.5*z0;
    double b=1.0/3.0+5.0/6.0*z0;
    g=min(a,b);
  
    }else if(z0>0.0 && z0<=0.5){

    double a=1100.0;
    double w=1.0/(1.0+a*square(square(z0-0.5)))/(1.0+a*square(square(z0-0.5)));
    double gp=(1.0/3.0+5.0/6.0*z0)*w+2.0*z0*(1.0-w);
    g=min(gp,2.0*z0);

    }else if(z0>0.5 && z0<1.0){

    double a=800.0;
    double w=1.0/(1.0+a*square(square(z0-0.5)))/(1.0+a*square(square(z0-0.5)));
    double gp=(1.0/3.0+5.0/6.0*z0)*w+(0.15*z0+0.85)*(1.0-w);
    g=min(gp,0.15*z0+0.85);

    }else{
    double a=1.2*(z0-1.0)+1.0;
    double b=1.0/3.0+5.0/6.0*z0;
    g=min(a,b);
    
    }
  
  }


  if(ROUND_Method==3){
  
  double a1=1.0+12.0*z0*z0;
  double a2=1.0+5.0*(z0-1.0)*(z0-1.0);
  double pl=1100.0*(z0-0.05)*(z0-0.05)*(z0-0.05)*(0.47-z0)*(0.47-z0)*(0.47-z0);
  double pr=18000.0*(z0-0.55)*(z0-0.55)*(z0-0.55)*(0.97-z0)*(0.97-z0)*(0.97-z0)*(0.97-z0)*(0.97-z0);
  double p1=5.0/6.0*z0+1.0/3.0+max(pl,0.0)+max(pr,0.0);
  double p2=1.5*z0;
  double p3=0.5*z0+0.5;
  double wc1=1.0/a1/a1/a1/a1;
  double wc2=1.0/a2/a2/a2/a2/a2/a2/a2/a2;

  g=(p1*(1.0-wc1)+p2*wc1)*(1.0-wc2)+p3*wc2;
  }


  return g*(u3-u1)+u1;

}

void reconst(double *rho_L, double *rho_R, double *rho, double *rhou_L, double *rhou_R, double *rhou, double *E_L, double *E_R, double *E)
{
	int i,j;
	double Lin_rho_R[nn],Lin_rho_L[nn];              // linear upwind
	double Lin_rhou_R[nn],Lin_rhou_L[nn];            // linear upwind
	double Lin_E_R[nn],Lin_E_L[nn];                  // linear upwind
	double ThiS_rho_R[nn],ThiS_rho_L[nn];            // THINC ��:small
	double ThiS_rhou_R[nn],ThiS_rhou_L[nn];          // THINC ��:small
	double ThiS_E_R[nn],ThiS_E_L[nn];                // THINC ��:small
	double ThiB_rho_R[nn],ThiB_rho_L[nn];            // THINC ��:large
	double ThiB_rhou_R[nn],ThiB_rhou_L[nn];          // THINC ��:large
	double ThiB_E_R[nn],ThiB_E_L[nn];                // THINC ��:large
	double WENO_rho_R[nn],WENO_rho_L[nn];            // WENO
	double WENO_rhou_R[nn],WENO_rhou_L[nn];          // WENO
	double WENO_E_R[nn],WENO_E_L[nn];  
	double ROUND_rho_R[nn],ROUND_rho_L[nn];            // WENO
	double ROUND_rhou_R[nn],ROUND_rhou_L[nn];          // WENO
	double ROUND_E_R[nn],ROUND_E_L[nn];                // WENO
	double US_rho_Rf[nn],US_rho_Lf[nn];              // select 1st
	double US_rhou_Rf[nn],US_rhou_Lf[nn];            // select 1st
	double US_E_Rf[nn],US_E_Lf[nn];                  // select 1st
	double US_rho_Rff[nn],US_rho_Lff[nn];            // select 2nd
	double US_rhou_Rff[nn],US_rhou_Lff[nn];          // select 2nd
	double US_E_Rff[nn],US_E_Lff[nn];                // select 2nd
	double w1_l[3],w2_l[3],w3_l[3],w4_l[3],w5_l[3];  // characteristic variable (left)
	double w1_r[3],w2_r[3],w3_r[3],w4_r[3],w5_r[3];  // characteristic variable (right)
	double Lin_w1_L[nn],Lin_w2_L[nn],Lin_w3_L[nn];
	double Lin_w1_R[nn],Lin_w2_R[nn],Lin_w3_R[nn];
	double WENO_w1_L[nn],WENO_w2_L[nn],WENO_w3_L[nn];
	double WENO_w1_R[nn],WENO_w2_R[nn],WENO_w3_R[nn];
	double ROUND_w1_L[nn],ROUND_w2_L[nn],ROUND_w3_L[nn];
	double ROUND_w1_R[nn],ROUND_w2_R[nn],ROUND_w3_R[nn];

	double u[nn],c[nn],H[nn];
	double u_ave[nn],c_ave[nn],H_ave[nn];

		double up[5];
	double un[5];
	
	set_u(u,c,H,rho,rhou,E);
	for(i=NGST-2;i!=NX+NGST+2;i++){
	/* Right */
		u_ave[i] = (sqrt(rho[i])*u[i]+sqrt(rho[i+1])*u[i+1]) / (sqrt(rho[i])+sqrt(rho[i+1]));
		H_ave[i] = (sqrt(rho[i])*H[i]+sqrt(rho[i+1])*H[i+1]) / (sqrt(rho[i])+sqrt(rho[i+1]));
		c_ave[i] = sqrt((gamma-1)*(H_ave[i]-u_ave[i]*u_ave[i]/2));

		eigen_L(u_ave[i],c_ave[i],rho[i-2],rhou[i-2],E[i-2],w1_r);
		eigen_L(u_ave[i],c_ave[i],rho[i-1],rhou[i-1],E[i-1],w2_r);
		eigen_L(u_ave[i],c_ave[i],rho[i],  rhou[i],  E[i],  w3_r);
		eigen_L(u_ave[i],c_ave[i],rho[i+1],rhou[i+1],E[i+1],w4_r);
		eigen_L(u_ave[i],c_ave[i],rho[i+2],rhou[i+2],E[i+2],w5_r);

		if(scheme == 1){
			WENO_w1_R[i]   = WENO(w1_r[0],w2_r[0],w3_r[0],w4_r[0],w5_r[0]);
			WENO_w2_R[i]   = WENO(w1_r[1],w2_r[1],w3_r[1],w4_r[1],w5_r[1]);
			WENO_w3_R[i]   = WENO(w1_r[2],w2_r[2],w3_r[2],w4_r[2],w5_r[2]);
		}
		else{
			ROUND_w1_R[i]   = ROUND(w2_r[0],w3_r[0],w4_r[0]);
			ROUND_w2_R[i]   = ROUND(w2_r[1],w3_r[1],w4_r[1]);
			ROUND_w3_R[i]   = ROUND(w2_r[2],w3_r[2],w4_r[2]);
		}

	/* Left */
		u_ave[i] = (sqrt(rho[i])*u[i]+sqrt(rho[i-1])*u[i-1]) / (sqrt(rho[i])+sqrt(rho[i-1]));
		H_ave[i] = (sqrt(rho[i])*H[i]+sqrt(rho[i-1])*H[i-1]) / (sqrt(rho[i])+sqrt(rho[i-1]));
		c_ave[i] = sqrt((gamma-1)*(H_ave[i]-u_ave[i]*u_ave[i]/2));

		eigen_L(u_ave[i],c_ave[i],rho[i+2],rhou[i+2],E[i+2],w1_l);
		eigen_L(u_ave[i],c_ave[i],rho[i+1],rhou[i+1],E[i+1],w2_l);
		eigen_L(u_ave[i],c_ave[i],rho[i],  rhou[i],  E[i],  w3_l);
		eigen_L(u_ave[i],c_ave[i],rho[i-1],rhou[i-1],E[i-1],w4_l);
		eigen_L(u_ave[i],c_ave[i],rho[i-2],rhou[i-2],E[i-2],w5_l);

		if(scheme == 1){
			WENO_w1_L[i] = WENO(w1_l[0],w2_l[0],w3_l[0],w4_l[0],w5_l[0]);
			WENO_w2_L[i] = WENO(w1_l[1],w2_l[1],w3_l[1],w4_l[1],w5_l[1]);
			WENO_w3_L[i] = WENO(w1_l[2],w2_l[2],w3_l[2],w4_l[2],w5_l[2]);
		}
		else{
			ROUND_w1_L[i] = ROUND(w2_l[0],w3_l[0],w4_l[0]);
			ROUND_w2_L[i] = ROUND(w2_l[1],w3_l[1],w4_l[1]);
			ROUND_w3_L[i] = ROUND(w2_l[2],w3_l[2],w4_l[2]);
		}

	}

		
	
	/*  convert to conservative */
	if(scheme==1){ // WENO
		for(i=NGST-1;i!=NX+NGST+1;i++){
		/* Right */
			u_ave[i] = (sqrt(rho[i])*u[i]+sqrt(rho[i+1])*u[i+1]) / (sqrt(rho[i])+sqrt(rho[i+1]));
			H_ave[i] = (sqrt(rho[i])*H[i]+sqrt(rho[i+1])*H[i+1]) / (sqrt(rho[i])+sqrt(rho[i+1]));
			c_ave[i] = sqrt((gamma-1)*(H_ave[i]-u_ave[i]*u_ave[i]/2));
			eigen_R(u_ave[i],c_ave[i],H_ave[i],&WENO_rho_R[i],&WENO_rhou_R[i],&WENO_E_R[i],WENO_w1_R[i],WENO_w2_R[i],WENO_w3_R[i]);
		/* Left */
			u_ave[i] = (sqrt(rho[i])*u[i]+sqrt(rho[i-1])*u[i-1]) / (sqrt(rho[i])+sqrt(rho[i-1]));
			H_ave[i] = (sqrt(rho[i])*H[i]+sqrt(rho[i-1])*H[i-1]) / (sqrt(rho[i])+sqrt(rho[i-1]));
			c_ave[i] = sqrt((gamma-1)*(H_ave[i]-u_ave[i]*u_ave[i]/2));
			eigen_R(u_ave[i],c_ave[i],H_ave[i],&WENO_rho_L[i],&WENO_rhou_L[i],&WENO_E_L[i],WENO_w1_L[i],WENO_w2_L[i],WENO_w3_L[i]);
		}
	}

	if(scheme==2){ // WENO3
		for(i=NGST-1;i!=NX+NGST+1;i++){
		/* Right */
			u_ave[i] = (sqrt(rho[i])*u[i]+sqrt(rho[i+1])*u[i+1]) / (sqrt(rho[i])+sqrt(rho[i+1]));
			H_ave[i] = (sqrt(rho[i])*H[i]+sqrt(rho[i+1])*H[i+1]) / (sqrt(rho[i])+sqrt(rho[i+1]));
			c_ave[i] = sqrt((gamma-1)*(H_ave[i]-u_ave[i]*u_ave[i]/2));
			eigen_R(u_ave[i],c_ave[i],H_ave[i],&ROUND_rho_R[i],&ROUND_rhou_R[i],&ROUND_E_R[i],ROUND_w1_R[i],ROUND_w2_R[i],ROUND_w3_R[i]);
		/* Left */
			u_ave[i] = (sqrt(rho[i])*u[i]+sqrt(rho[i-1])*u[i-1]) / (sqrt(rho[i])+sqrt(rho[i-1]));
			H_ave[i] = (sqrt(rho[i])*H[i]+sqrt(rho[i-1])*H[i-1]) / (sqrt(rho[i])+sqrt(rho[i-1]));
			c_ave[i] = sqrt((gamma-1)*(H_ave[i]-u_ave[i]*u_ave[i]/2));
			eigen_R(u_ave[i],c_ave[i],H_ave[i],&ROUND_rho_L[i],&ROUND_rhou_L[i],&ROUND_E_L[i],ROUND_w1_L[i],ROUND_w2_L[i],ROUND_w3_L[i]);
		}
	}
	
	if(scheme==1){
		for(i=NGST;i!=NX+NGST+1;i++){
			rho_L[i]  = WENO_rho_R[i-1];
			rho_R[i]  = WENO_rho_L[i];
			rhou_L[i] = WENO_rhou_R[i-1];
			rhou_R[i] = WENO_rhou_L[i];
			E_L[i]    = WENO_E_R[i-1];
			E_R[i]    = WENO_E_L[i];
		}
	}else{
		for(i=NGST;i!=NX+NGST+1;i++){
			rho_L[i]  = ROUND_rho_R[i-1];
			rho_R[i]  = ROUND_rho_L[i];
			rhou_L[i] = ROUND_rhou_R[i-1];
			rhou_R[i] = ROUND_rhou_L[i];
			E_L[i]    = ROUND_E_R[i-1];
			E_R[i]    = ROUND_E_L[i];
		}
	}
}

double square(double a)
{
  return a*a;
}
