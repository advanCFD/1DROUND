#include "define.h"
#include "function.h"

void EF(int n, double *dt, double dx, double CFL, double *rho, double *rhou, double *E, double *rho_b, double *rhou_b, double *E_b, double t, double total_t)

{
	int i;
	double rho_s[nn], rhou_s[nn], E_s[nn];
	double rho_l[nn],rhou_l[nn],E_l[nn];                 // �Z�����E�l(��)
	double rho_r[nn],rhou_r[nn],E_r[nn];                 // �Z�����E�l(�E)
	double rho_Lf[nn],rhou_Lf[nn],E_Lf[nn];              // ���U���I�y���[�^
	
	for(i=0;i!=nn;i++){
		rho_Lf[i] = 0.0;
		rhou_Lf[i] = 0.0;
		E_Lf[i] = 0.0;
	}
	if(bc == 2){
		wall(rho_b,rhou_b,E_b);
	}else{
		neumann(rho_b);
		neumann(rhou_b);
		neumann(E_b);
	}
	reconst(rho_l,rho_r,rho_b,rhou_l,rhou_r,rhou_b,E_l,E_r,E_b);
	if(riemann_sol == 1) riemann_hll(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
	else if(riemann_sol == 2) riemann_hllc(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
	discre_L(dx,*dt,rho_s,rho_Lf);
	discre_L(dx,*dt,rhou_s,rhou_Lf);
	discre_L(dx,*dt,E_s,E_Lf);
	
	for(i=NGST;i!=NX+NGST;i++){
		rho[i] = rho_b[i] + rho_Lf[i];
		rhou[i] = rhou_b[i] + rhou_Lf[i];
		E[i] = E_b[i] + E_Lf[i];
	}
}


void SSP_RK3(int n, double *dt, double dx, double CFL, double *rho, double *rhou, double *E, double *rho_b, double *rhou_b, double *E_b, double t, double total_t)
{
	int i, rl;
	double rho_s[nn], rhou_s[nn], E_s[nn];
	double rho_l[nn],u_l[nn],p_l[nn],rhou_l[nn],E_l[nn];    // �Z�����E�l(��)
	double rho_r[nn],u_r[nn],p_r[nn],rhou_r[nn],E_r[nn];    // �Z�����E�l(�E)
	double rho_Lf[nn],rhou_Lf[nn],E_Lf[nn];                 // ���U���I�y���[�^

	for(rl = 0; rl < 3 ; rl++){
		if(bc == 2){
			if(rl == 0) wall(rho_b,rhou_b,E_b);
			else wall(rho,rhou,E);
		}else{
			if(rl == 0){
				neumann(rho_b);
				neumann(rhou_b);
				neumann(E_b);
			}else{
				neumann(rho);
				neumann(rhou);
				neumann(E);
			}
		}
		if(rl==0) reconst(rho_l,rho_r,rho_b,rhou_l,rhou_r,rhou_b,E_l,E_r,E_b);
		else      reconst(rho_l,rho_r,rho,  rhou_l,rhou_r,rhou,  E_l,E_r,E);
		if(riemann_sol == 1)      riemann_hll(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
		else if(riemann_sol == 2) riemann_hllc(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
		discre_L(dx,*dt,rho_s,rho_Lf);
		discre_L(dx,*dt,rhou_s,rhou_Lf);
		discre_L(dx,*dt,E_s,E_Lf);
		if(rl==0){
			for(i=NGST;i!=NX+NGST;i++){
				rho[i] = rho_b[i] + rho_Lf[i];
				rhou[i] = rhou_b[i] + rhou_Lf[i];
				E[i] = E_b[i] + E_Lf[i];
			}
		}else{
			for(i=NGST;i!=NX+NGST;i++){
				rho[i] = rho[i] + rho_Lf[i];
				rhou[i] = rhou[i] + rhou_Lf[i];
				E[i] = E[i] + E_Lf[i];
			}
		}
		if(rl==1){
			for(i=NGST;i!=NX+NGST;i++){
				rho[i] = 3.0 * rho_b[i] / 4.0 + 1.0 * rho[i] / 4.0;
				rhou[i] = 3.0 * rhou_b[i] / 4.0 + 1.0 * rhou[i] / 4.0;
				E[i] = 3.0 * E_b[i] / 4.0 + 1.0 * E[i] / 4.0;
			}
		}
		if(rl==2){
			for(i=NGST;i!=NX+NGST;i++){
				rho[i] = 1.0 * rho_b[i] / 3.0 + 2.0 * rho[i] / 3.0;
				rhou[i] = 1.0 * rhou_b[i] / 3.0 + 2.0 * rhou[i] / 3.0;
				E[i] = 1.0 * E_b[i] / 3.0 + 2.0 * E[i] / 3.0;
			}
		}
	}
}



void RK5(int n, double *dt, double dx, double CFL, double *rho, double *rhou, double *E, double *rho_b, double *rhou_b, double *E_b, double t, double total_t)
{
	int i,rl;
	double rho_b1[nn], rhou_b1[nn], E_b1[nn];
	double rho_b2[nn], rhou_b2[nn], E_b2[nn];
	double rho_b3[nn], rhou_b3[nn], E_b3[nn];
	double rho_b4[nn], rhou_b4[nn], E_b4[nn];
	double rho_s[nn], rhou_s[nn], E_s[nn];
	double rho_l[nn],u_l[nn],p_l[nn],rhou_l[nn],E_l[nn];         // �Z�����E�l(��)
	double rho_r[nn],u_r[nn],p_r[nn],rhou_r[nn],E_r[nn];         // �Z�����E�l(�E)
	double rho_Lf[nn],rhou_Lf[nn],E_Lf[nn];                      // ���U���I�y���[�^
	
	for(rl = 0; rl < 5 ; rl++){
		if(rl==0){
			if(bc == 2){
				wall(rho_b,rhou_b,E_b);
			}else{
				neumann(rho_b);
				neumann(rhou_b);
				neumann(E_b);
			}
			reconst(rho_l,rho_r,rho_b,rhou_l,rhou_r,rhou_b,E_l,E_r,E_b);
			if(riemann_sol == 1) riemann_hll(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
			else if(riemann_sol == 2) riemann_hllc(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
			discre_L(dx,*dt,rho_s,rho_Lf);
			discre_L(dx,*dt,rhou_s,rhou_Lf);
			discre_L(dx,*dt,E_s,E_Lf);
			for(i=NGST;i!=NX+NGST;i++){
				rho_b1[i] = rho_b[i]+0.391752226571890*rho_Lf[i];
				rhou_b1[i] = rhou_b[i]+0.391752226571890*rhou_Lf[i];
				E_b1[i] = E_b[i]+0.391752226571890*E_Lf[i];
			}
		}
		if(rl==1){
			if(bc == 2){
				wall(rho_b1,rhou_b1,E_b1);
			}else{
				neumann(rho_b1);
				neumann(rhou_b1);
				neumann(E_b1);
			}
			reconst(rho_l,rho_r,rho_b1,rhou_l,rhou_r,rhou_b1,E_l,E_r,E_b1);
			if(riemann_sol == 1) riemann_hll(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
			else if(riemann_sol == 2) riemann_hllc(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
			discre_L(dx,*dt,rho_s,rho_Lf);
			discre_L(dx,*dt,rhou_s,rhou_Lf);
			discre_L(dx,*dt,E_s,E_Lf);
			for(i=NGST;i!=NX+NGST;++i){
				rho_b2[i]=0.444370493651220*rho_b[i]+0.555629506348765*rho_b1[i]+0.368410593050371*rho_Lf[i];
				rhou_b2[i]=0.444370493651220*rhou_b[i]+0.555629506348765*rhou_b1[i]+0.368410593050371*rhou_Lf[i];
				E_b2[i]=0.444370493651220*E_b[i]+0.555629506348765*E_b1[i]+0.368410593050371*E_Lf[i];
			}
		}
		if(rl==2){
			if(bc == 2){
				wall(rho_b2,rhou_b2,E_b2);
			}else{
				neumann(rho_b2);
				neumann(rhou_b2);
				neumann(E_b2);
			}
			reconst(rho_l,rho_r,rho_b2,rhou_l,rhou_r,rhou_b2,E_l,E_r,E_b2);
			if(riemann_sol == 1) riemann_hll(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
			else if(riemann_sol == 2) riemann_hllc(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
			discre_L(dx,*dt,rho_s,rho_Lf);
			discre_L(dx,*dt,rhou_s,rhou_Lf);
			discre_L(dx,*dt,E_s,E_Lf);
			for(i=NGST;i!=NX+NGST;++i){
				rho_b3[i]=0.620101851488403*rho_b[i]+0.379898148511597*rho_b2[i]+0.251891774271694*rho_Lf[i];
				rhou_b3[i]=0.620101851488403*rhou_b[i]+0.379898148511597*rhou_b2[i]+0.251891774271694*rhou_Lf[i];
				E_b3[i]=0.620101851488403*E_b[i]+0.379898148511597*E_b2[i]+0.251891774271694*E_Lf[i];
			}
		}
		if(rl==3){
			if(bc == 2){
				wall(rho_b3,rhou_b3,E_b3);
			}else{
				neumann(rho_b3);
				neumann(rhou_b3);
				neumann(E_b3);
			}
			reconst(rho_l,rho_r,rho_b3,rhou_l,rhou_r,rhou_b3,E_l,E_r,E_b3);
			if(riemann_sol == 1) riemann_hll(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
			else if(riemann_sol == 2) riemann_hllc(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
			discre_L(dx,*dt,rho_s,rho_Lf);
			discre_L(dx,*dt,rhou_s,rhou_Lf);
			discre_L(dx,*dt,E_s,E_Lf);
			for(i=NGST;i!=NX+NGST;++i){
				rho_b4[i]=0.178079954393132*rho_b[i]+0.821920045606868*rho_b3[i]+0.544974750228521*rho_Lf[i];
				rhou_b4[i]=0.178079954393132*rhou_b[i]+0.821920045606868*rhou_b3[i]+0.544974750228521*rhou_Lf[i];
				E_b4[i]=0.178079954393132*E_b[i]+0.821920045606868*E_b3[i]+0.544974750228521*E_Lf[i];
			}
		}
		if(rl==4){
			if(bc == 2){
				wall(rho_b4,rhou_b4,E_b4);
			}else{
				neumann(rho_b4);
				neumann(rhou_b4);
				neumann(E_b4);
			}
			reconst(rho_l,rho_r,rho_b4,rhou_l,rhou_r,rhou_b4,E_l,E_r,E_b4);
			if(riemann_sol == 1) riemann_hll(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
			else if(riemann_sol == 2) riemann_hllc(dt,dx,CFL,rho_l,rho_r,rho_s,rhou_l,rhou_r,rhou_s,E_l,E_r,E_s,t,total_t);
			discre_L(dx,*dt,rho_s,rho_Lf);
			discre_L(dx,*dt,rhou_s,rhou_Lf);
			discre_L(dx,*dt,E_s,E_Lf);
			for(i=NGST;i!=NX+NGST;++i){
				rho[i]=0.517231671970585*rho_b2[i]-0.020812619136066*rho_b[i]+0.503580947165482*rho_b4[i]+0.226007483236906*rho_Lf[i];
				rhou[i]=0.517231671970585*rhou_b2[i]-0.020812619136066*rhou_b[i]+0.503580947165482*rhou_b4[i]+0.226007483236906*rhou_Lf[i];
				E[i]=0.517231671970585*E_b2[i]-0.020812619136066*E_b[i]+0.503580947165482*E_b4[i]+0.226007483236906*E_Lf[i];
			}
		}
	}
}

