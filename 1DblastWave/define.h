#define GNUPLOT_PATH "C:/PROGRA~2/gnuplot/bin/gnuplot.exe -persist"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <iostream>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define max3(a,b,c) ((a)>max(b,c)?(a):(max(b,c)))
#define min3(a,b,c) ((a)<min(b,c)?(a):(min(b,c)))
#define pi M_PI

#define scheme 2       // 1:WENO5 scheme 2:ROUND scheme
#define WENO_Method 0  // 0:WENO-JS 1:WENO-Z
#define ROUND_Method 0  // 0:ROUND_A+ 1:ROUND_W 2:ROUND_L  3:ROUND_L+  
#define time_in 3      // 1:euler 2:RK3 3:RK5,4
#define riemann_sol 2  // 1:HLL 2:HLLC
#define bc 2           // 1:neumann 2:wall

#define ic 3           // 1: sod   2:lax   3: blast wave  4: density interaction 5: Le Blanc 6:Sedov 
#define NX 400         // sod,lax 100;    blast wave, density interaction 400;     Leblanc 200  Sedov 801; 
// 600 cells for density interaction, WENO3, final_t=0.18

#define beta_S 1.1     // Value of small
#define beta_B 1.6     // Value of large 
#define method 1       // 1:thinc 2:TVD
#define lim 4          // 1:UW 2:minmod 3:superbee 4:VanLeer

#define NGST 5         // Ghost cell
#define nn NX+2*NGST   // All cell
#define xlow 0.0     //  Sedov [-2.0025:2.0025】 
#define xup 1.0
#define eps 1e-6
#define epsilon 1e-6
#define epsilon1 1e-20

#define gamma 1.4


