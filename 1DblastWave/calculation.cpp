#include "define.h"
#include "function.h"


void cal_up(double *rho, double *rhou, double *E, double *u, double *p)
{
	int i;
	
	for(i=NGST-2;i<nn;i++){
		u[i] = rhou[i] / rho[i];
		p[i] = (gamma-1.0) * (E[i] - 0.5*u[i]*u[i]*rho[i]);
	}
}


void cal_rhouE(double *rho, double *rhou, double *E, double *u, double *p)
{
	int i;

	for(i=NGST-2;i<nn;i++){
		rhou[i] = u[i] * rho[i];
		E[i] = p[i] / (gamma-1.0) + 0.5 * rho[i] * u[i]*u[i];
	}
}


void set_u(double *u, double *c, double *H, double *rho, double *rhou, double *E)
{
	int i;
	double p;
	
	for(i=0;i<nn;i++){
		u[i] = rhou[i] / rho[i];
		p = (gamma-1.0) * (E[i] - 0.5*u[i]*u[i]*rho[i]);
		c[i] = sqrt(gamma * p / rho[i]);
		H[i] = (E[i]+p)/rho[i];
	}
}



void eigen_L(double u, double c, double rho, double rhou, double E, double *w)
{
	double a1,a2;
	
	a1 = u*u/2*(gamma-1)/(c*c);
	a2 = (gamma-1)/(c*c);
	w[0] = (a1+u/c)*rho/2 - (1/c+a2*u)*rhou/2 + a2*E/2;
	w[1] = (1-a1)*rho     + a2*u*rhou         - a2*E;
	w[2] = (a1-u/c)*rho/2 + (1/c-a2*u)*rhou/2 + a2*E/2;
}


void eigen_R(double u, double c, double H, double *rho, double *rhou, double *E, double w1, double w2, double w3)
{
	
	*rho  = w1         + w2       + w3;
	*rhou = (u-c)*w1   + u*w2     + (u+c)*w3;
	*E    = (H-u*c)*w1 + u*u*w2/2 + (H+u*c)*w3;
}


void discre_L(double dx,double dt, double *flux, double *Lf)
{
	int i;

	for(i=NGST;i<=NX+NGST;i++) Lf[i] = dt / dx * (flux[i] - flux[i+1]);
}


void sound_speed(double rho, double rhou, double E, double *p, double *c)
{
	double e,u;
	
	u = rhou / rho;
	e = E - 0.5*u*u*rho;
	*p = (gamma-1.0) * e;
	*c = sqrt(gamma * (*p) / rho);
}


void calflux(double rho, double rhou, double E, double p, double u, double *rho_flux, double *rhou_flux, double *E_flux)
{
	*rho_flux = rho * u;
	*rhou_flux = *rho_flux * u + p;
	*E_flux = u * (E + p);
}


void cal_dt(double *s_l, double *s_r, double *dt, double dx, double CFL , double t, double total_t)
{
	int i;
	double umax,u,cfl_local;

	umax = 0;
	for(i=NGST-1;i<NX+NGST;i++){
		u = max(fabs(s_l[i]),fabs(s_r[i]));
		if(u >= umax) umax = u;
	}
	*dt = CFL*dx/umax;
	if((t+(*dt) > total_t) && (t != total_t)){
      		*dt = total_t - t;
      		
  	 }
}


void m_state(double *rhom, double *rhoum, double *Em, double s_m, double s, double rho_flux, double rhou_flux, double E_flux, double rho, double rhou, double E, double p)
{
	double rho_k,rhou_k,E_k,pm,uk;
   
	uk = rhou / rho;
	rho_k = s * rho - rho_flux;
	rhou_k = s * rhou - rhou_flux;
	E_k = s * E - E_flux;
	*rhom = rho_k / (s - s_m);
	*rhoum = (*rhom) * s_m;
	pm = s_m * rho_k - rhou_k;
	*Em = (E_k + s_m * pm) / (s-s_m);  
}


double HLL(double s_l, double s_r, double flux_l, double flux_r, double f_l, double f_r)
{
	return  ( (s_r*flux_l - s_l*flux_r + s_l*s_r*(f_r-f_l)) / (s_r-s_l) );
}


void riemann_fan(double *f_fan, double s_l, double s_r, double f_l, double f_r, double fl_flux, double fr_flux)
{   
	*f_fan = ((s_r*f_r - s_l*f_l) - (fr_flux - fl_flux)) / (s_r-s_l);
}


double HLLC(double s, double flux, double f, double fm)
{
      return ( flux + s * (fm-f) );
}
