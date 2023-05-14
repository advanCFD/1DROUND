
/* riemann */
void riemann_hll(double *dt, double dx, double CFL, double *rho_ls, double *rho_rs, double *rho_s, double *rhou_ls, double *rhou_rs, double *rhou_s, double *E_ls, double *E_rs, double *E_s, double t, double total_t);
void riemann_hllc(double *dt, double dx, double CFL, double *rho_ls, double *rho_rs, double *rho_s, double *rhou_ls, double *rhou_rs, double *rhou_s, double *E_ls, double *E_rs, double *E_s, double t, double total_t);

/* time_integral */
void EF(int n, double *dt, double dx, double CFL, double *rho, double *rhou, double *E, double *rho_b, double *rhou_b, double *E_b, double t, double total_t);
void SSP_RK3(int n, double *dt, double dx, double CFL, double *rho, double *rhou, double *E, double *rho_b, double *rhou_b, double *E_b, double t, double total_t);
void RK5(int n, double *dt, double dx, double CFL, double *rho, double *rhou, double *E, double *rho_b, double *rhou_b, double *E_b, double t, double total_t);

/* reconstruction */
void select_function_1st(double *ThiM_L,double *ThiM_R, double *Lin_L, double *Lin_R, double *US_Lf, double *US_Rf);
void select_function_2nd(double *ThiS_L,double *ThiS_R, double *US_Lf, double *US_Rf, double *US_Lff, double *US_Rff);
void select_function_3rd(double *ThiB_L, double *ThiB_R, double *US_Lff, double *US_Rff, double *US_Lfff, double *US_Rfff);
double Linear(double w1, double w2, double w3, double w4, double w5, double w6, double w7);
void thinc(int direction, double a, double b, double c, double *fs, double beta);
double WENO(double w1, double w2, double w3, double w4, double w5);
double ROUND(double w1, double w2, double w3);
void reconst(double *rho_L, double *rho_R, double *rho, double *rhou_L, double *rhou_R, double *rhou, double *E_L, double *E_R, double *E);

/* boundary_condition */
void period(double *f);
void wall(double *rho, double *rhou, double *E);
void neumann(double *f);


/* calculation */
void cal_up(double *rho, double *rhou, double *E, double *u, double *p);
void cal_rhouE(double *rho, double *rhou, double *E, double *u, double *p);
void set_u(double *u_ave, double *c_ave, double *H_ave, double *rho, double *rhou, double *E);
void eigen_L(double u, double c, double rho, double rhou, double E, double *w);
void eigen_R(double u, double c, double H, double *rho, double *rhou, double *E, double w1, double w2, double w3);
void discre_L(double dx, double dt, double *flux, double *Lf);
void sound_speed(double rho, double rhou, double E, double *p, double *c);
void calflux(double rho, double rhou, double E, double p, double u, double *rho_flux, double *rhou_flux, double *E_flux);
void cal_dt(double *s_l, double *s_r, double *dt, double dx, double CFL, double t, double total_t);
void m_state(double *rhom, double *rhoum, double *Em, double s_m, double s, double rho_flux, double rhou_flux, double E_flux, double rho, double rhou, double E, double p);
double HLL(double s_l, double s_r, double flux_l, double flux_r, double f_l, double f_r);
void riemann_fan(double *f_fan, double s_l, double s_r, double f_l, double f_r, double fl_flux, double fr_flux);
double HLLC(double s, double flux, double f, double fm);
double square(double a);

/* initial_condition */
void initial(double *x, double dx, double *rho, double *rhou, double *E, double *rho_b, double *rhou_b, double *E_b);

/* write */
void write(FILE *gp, int n, double *x, double *rho, double *rhou, double *E, int output_step);

/* plot */
void make_anime(FILE *fp);
void read_exact(FILE *ep, double *xe, double *fe, double *pe, double *ue);
void graph_exact(FILE *fp, double *x, double *fc, double *xe, double *fe);
void graph(FILE *fp, double *x, double *fc);
void Output(double *x, double *rho, double *rhou, double *E, int step);

void thinc2(int direction, double a, double b, double c, double *fs, double beta);
void thinc_try(int direction, double a, double b, double c, double d, double e, double *fs);
