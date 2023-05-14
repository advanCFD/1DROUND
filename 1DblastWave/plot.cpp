#include "define.h"
#include "function.h"

void make_anime(FILE *fp)
{
	fprintf(fp,"set terminal windows\n");
	fprintf(fp,"set terminal gif animate optimize size \n");
    fprintf(fp,"set output 'anime.gif'\n");
}

void read_exact(FILE *ep, double *xe, double *fe, double *pe, double *ue)
{
	int i,data;
	char fname[25];
	
	data = 1000;
	sprintf(fname,"exact_shock.csv");
	
	ep = fopen(fname,"r");
	if( ep == NULL ){
    		printf( "Can't open %s\n", fname);
    		exit(0);
  	}
	for(i=0;i<data;i++){
		fscanf(ep, "%lf,%lf,%lf,%lf\n", &xe[i],&fe[i],&pe[i],&ue[i]);
    		//printf("%d %lf %lf\n",i,xe[i],fe[i]);
  	}
	fclose(ep);
}


void graph_exact(FILE *fp, double *x, double *fc, double *xe, double *fe)
{
	int i,ret,data;
	
	fprintf(fp, "set multiplot\n");
	fprintf(fp, "unset key \n");
	fprintf(fp, "set format y \"%%f\" \n");
	fprintf(fp, "set xrange [-5.0:5.0] \n");
	
	data = 1000;
	fprintf(fp, "set yrange [0.8:1.7] \n");
	
	fprintf(fp, "plot '-' with l lt -1 lw 2 \n");
	for(i=0;i<data;i++) fprintf(fp,"%f\t%f\t\n",xe[i],fe[i]);
	fprintf(fp, "e \n");

	fprintf(fp, "plot (10e-10*x), '-' with lp lt 7 lw 2 pt 7 ps 1 lc 1\n");
	for(i=NGST;i<NX+NGST;i++) fprintf(fp,"%f\t%f\t\n",x[i],fc[i]);
	fprintf(fp, "e \n");

	fprintf(fp, "set nomultiplot\n");
	fflush(fp);
}


void graph(FILE *fp, double *x, double *fc)
{
	int i;
	fprintf(fp, "unset key \n");
	fprintf(fp, "set format y \"%%f\" \n");
	fprintf(fp, "set yrange [0.8:1.7] \n");
	fprintf(fp, "set xrange [-5.0:5.0] \n");

	fprintf(fp, "plot (10e-10*x), '-' with lp lt 7 lw 2 pt 7 ps 1 lc 1 \n");
	for(i=NGST;i<NX+NGST;i++) fprintf(fp,"%f\t%f\t\n",x[i],fc[i]);
	fprintf(fp, "e \n");

	fflush(fp);
}
