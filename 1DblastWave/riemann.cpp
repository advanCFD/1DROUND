#include "define.h"
#include "function.h"


void riemann_hll(double *dt, double dx, double CFL, double *rho_ls, double *rho_rs, double *rho_s, double *rhou_ls, double *rhou_rs, double *rhou_s, double *E_ls, double *E_rs, double *E_s, double t, double total_t)
{
	int i;
	double u_l,u_r;
	double p_l,p_r;
	double c_l,c_r;
	double s_l[nn],s_r[nn];
	double rhol_flux, rhoul_flux, El_flux;
	double rhor_flux, rhour_flux, Er_flux;
	
	for(i=NGST-1;i<NX+NGST+1;i++){
		sound_speed(rho_ls[i],rhou_ls[i],E_ls[i],&p_l,&c_l);
		sound_speed(rho_rs[i],rhou_rs[i],E_rs[i],&p_r,&c_r);
		u_l = rhou_ls[i] / rho_ls[i];
		u_r = rhou_rs[i] / rho_rs[i];
		s_l[i] = min(u_l - c_l, u_r - c_r);
		s_r[i] = max(u_l + c_l, u_r + c_r);
		calflux(rho_ls[i],rhou_ls[i],E_ls[i],p_l,u_l,&rhol_flux,&rhoul_flux,&El_flux);
		calflux(rho_rs[i],rhou_rs[i],E_rs[i],p_r,u_r,&rhor_flux,&rhour_flux,&Er_flux);
		if(s_l[i]>=0.0){
			rho_s[i] = rhol_flux;
			rhou_s[i] = rhoul_flux;
			E_s[i] = El_flux;
		}else if(s_l[i]<=0.0 && s_r[i]>=0.0){
			rho_s[i] = HLL(s_l[i],s_r[i],rhol_flux,rhor_flux,rho_ls[i],rho_rs[i]);
			rhou_s[i] = HLL(s_l[i],s_r[i],rhoul_flux,rhour_flux,rhou_ls[i],rhou_rs[i]);
			E_s[i] = HLL(s_l[i],s_r[i],El_flux,Er_flux,E_ls[i],E_rs[i]);
		}else if(s_r[i]<=0.0){
			rho_s[i] = rhor_flux;
			rhou_s[i] = rhour_flux;
			E_s[i] = Er_flux;
		}
	}
	cal_dt(s_l,s_r,dt,dx,CFL,t,total_t);
}


void riemann_hllc(double *dt, double dx, double CFL, double *rho_ls, double *rho_rs, double *rho_s, double *rhou_ls, double *rhou_rs, double *rhou_s, double *E_ls, double *E_rs, double *E_s, double t, double total_t)
{
	int i;
	double u_l,u_r;
	double p_l,p_r;
	double c_l,c_r;
	double s_l[nn],s_r[nn],s_m[nn];
	double rhol_flux, rhoul_flux, El_flux;
	double rhor_flux, rhour_flux, Er_flux;
	double rho_fan, rhou_fan, E_fan;
	double rhom_l, rhoum_l, Em_l;
	double rhom_r, rhoum_r, Em_r;
	
	for(i=NGST;i<NX+NGST+1;i++){
		sound_speed(rho_ls[i],rhou_ls[i],E_ls[i],&p_l,&c_l);
		sound_speed(rho_rs[i],rhou_rs[i],E_rs[i],&p_r,&c_r);
		u_l = rhou_ls[i] / rho_ls[i];
		u_r = rhou_rs[i] / rho_rs[i];
		s_l[i] = min(u_l - c_l, u_r - c_r);
		s_r[i] = max(u_l + c_l, u_r + c_r);
		calflux(rho_ls[i],rhou_ls[i],E_ls[i],p_l,u_l,&rhol_flux,&rhoul_flux,&El_flux);
		calflux(rho_rs[i],rhou_rs[i],E_rs[i],p_r,u_r,&rhor_flux,&rhour_flux,&Er_flux);
		riemann_fan(&rho_fan,s_l[i],s_r[i],rho_ls[i],rho_rs[i],rhol_flux,rhor_flux);
		riemann_fan(&rhou_fan,s_l[i],s_r[i],rhou_ls[i],rhou_rs[i],rhoul_flux,rhour_flux);
		riemann_fan(&E_fan,s_l[i],s_r[i],E_ls[i],E_rs[i],El_flux,Er_flux);
		s_m[i] = rhou_fan / rho_fan;
		if(s_l[i]>=0.0){
			rho_s[i] = rhol_flux;
			rhou_s[i] = rhoul_flux;
			E_s[i] = El_flux;
		}else if(s_l[i]<=0.0 && s_m[i]>=0.0){
			m_state(&rhom_l,&rhoum_l,&Em_l,s_m[i],s_l[i],rhol_flux,rhoul_flux,El_flux,rho_ls[i],rhou_ls[i],E_ls[i],p_l);
			rho_s[i] = HLLC(s_l[i],rhol_flux,rho_ls[i],rhom_l);
			rhou_s[i] = HLLC(s_l[i],rhoul_flux,rhou_ls[i],rhoum_l);
			E_s[i] = HLLC(s_l[i],El_flux,E_ls[i],Em_l);
		}else if(s_m[i]<=0.0 && s_r[i]>=0.0){
			m_state(&rhom_r,&rhoum_r,&Em_r,s_m[i],s_r[i],rhor_flux,rhour_flux,Er_flux,rho_rs[i],rhou_rs[i],E_rs[i],p_r);
			rho_s[i] = HLLC(s_r[i],rhor_flux,rho_rs[i],rhom_r);
			rhou_s[i] = HLLC(s_r[i],rhour_flux,rhou_rs[i],rhoum_r);
			E_s[i] = HLLC(s_r[i],Er_flux,E_rs[i],Em_r);
		}else if(s_r[i]<=0.0){
			rho_s[i] = rhor_flux;
			rhou_s[i] = rhour_flux;
			E_s[i] = Er_flux;
		}
	}
	cal_dt(s_l,s_r,dt,dx,CFL,t,total_t);
}
