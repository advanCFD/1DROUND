// ROUND reconstruction function in a compact stencil:
// |---u1---|---u2---|---u3---|   u1: upwind cell, u2: target cell, u3 downwind cell
//                   ^
//return:            ul           ul: reconstructed value at the left-side cell boundary 
//Details see Deng, Xi. "A Unified Framework for Non-linear Reconstruction Schemes in a Compact Stencil. Part 1: Beyond Second Order." Journal of Computational Physics (2023): 112052.

//low-dissipation, structure-preserving ROUND scheme see Eq (7.7):
double ROUNDAplus(double u1, double u2, double u3) 
{

  double zeps  = 1.0e-20;

  double z0=(u2-u1+zeps)/(u3-u1+zeps);
  
  double g=z0;
  
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

  return g*(u3-u1)+u1;

}

//dissipative ROUND scheme with a continuous function Eq. 6.6. 
double ROUNDW(double u1, double u2, double u3) 
{

  double zeps  = 1.0e-20;

  double z0=(u2-u1+zeps)/(u3-u1+zeps);

  double w1=1.0/square(1.0+25.0*z0*z0);
  
  double w2=1.0/square(1.0+25.0*(1.0-z0)*(1.0-z0));
  
  double g=(5.0/6.0+2.0/3.0*w1-1.0/3.0*w2)*z0+1.0/3.0-1.0/3.0*w1+1.0/6.0*w2;
  
  return g*(u3-u1)+u1;

}

//low-dissipation ROUND: improve the capability of resolving critical points. see Eq (7.1) amd Eq (7.7)
double ROUNDL(double u1, double u2, double u3)  
{

  double zeps  = 1.0e-20;

  double z0=(u2-u1+zeps)/(u3-u1+zeps);
  
  double g=z0;
  
  
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

  return g*(u3-u1)+u1;

}


// To further reduce dissipation. See Eq (7.5)
double ROUNDLplus(double u1, double u2, double u3) 
{

  double zeps  = 1.0e-20;

  double z0=(u2-u1+zeps)/(u3-u1+zeps);
  
  double g=z0;
  
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

  return g*(u3-u1)+u1;

}

//alternative formulation to Eq.(7.5).    See Eq.(7.4)
/*double ROUNDLplus(double u1, double u2, double u3) 
{

  double zeps  = 1.0e-20;

  double z0=(u2-u1+zeps)/(u3-u1+zeps);
  
  double g=z0;
  
  if(z0<=0.0){
   double a1=1.0+5.0*z0*z0;
   double w1=1.0-1.0/a1/a1/a1/a1;

   g=(1.0/3.0+5.0/6.0*z0)*w1;

   }else if(z0>0.0 && z0<=0.5){

   double a1=1.0+300.0*(z0-0.5)*(z0-0.5)*(z0-0.5)*(z0-0.5);

   double w1=1.0/a1/a1;

   double p2=2.5*z0;
  
   double p=(1.0/3.0+5.0/6.0*z0)*w1+p2*(1.0-w1);

   g=min(p,p2);

   }else if(z0>0.5 && z0<=1.0){
   double a1=1.0+800.0*(z0-0.5)*(z0-0.5)*(z0-0.5)*(z0-0.5);
   double w1=1.0/a1/a1;

   double p2=0.05*(z0-1.0)+1.0;
  
   double p=(1.0/3.0+5.0/6.0*z0)*w1+p2*(1.0-w1);

   g=min(p,p2);


  }else{
   double a1=1.0+20.0*(z0-1.0)*(z0-1.0);
   double w1=1.0/a1/a1;

   g=(1.0/3.0+5.0/6.0*z0)*(1.0-w1)+z0*w1;

  }

  return g*(u3-u1)+u1;

}*/
