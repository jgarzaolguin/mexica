#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

 double numerical_int(double   *grid,
                      double   *grid_fun,
                      int       points)
 
 {//function

 int unsigned i;

//////////////////////

 double suma1;
 double rho_0;
 double rho_1;
 double rho_2;
 double x0;
 double x1;
 double x2;
 double a;
 double b;
 double c;
 double d;
 double pow2_x0;
 double pow2_x1;
 double pow2_x2;
 double coef13;
 double coef12;
 int iglobal;

 coef13   = 1.f/3.f;
 coef12   = 1.f/2.f;

 suma1=0.f;

      for (i = 0; i < (points/2 - 1); i++) {
         iglobal = 2*(i + 1);

             x0 = grid[iglobal -2];           
             x1 = grid[iglobal -1];           
             x2 = grid[iglobal];              

             rho_0 = grid_fun[iglobal -2];
             rho_1 = grid_fun[iglobal -1] ;
             rho_2 = grid_fun[iglobal];

      rho_0  = x0*x0*rho_0;     //f(x0)
      rho_1  = x1*x1*rho_1;     //f(x1)
      rho_2  = x2*x2*rho_2;     //f(x2)

      pow2_x0 = pow(x0,2.f);
      pow2_x1 = pow(x1,2.f);
      pow2_x2 = pow(x2,2.f);

      d = (x0 - x1)*(x0 - x2)*(x1 - x2);

      a = ((x1 - x2)*rho_0 + (x2 - x0)*rho_1 + (x0 - x1)*rho_2)/d;

      b = ((rho_1 - rho_0)*pow2_x0 + (rho_2 - rho_0)*pow2_x1 + (rho_0 - rho_1)*pow2_x2)/d;

      c = ((rho_2*x1 - rho_1*x2)*pow2_x0 + (rho_0*x2 - rho_2*x0)*pow2_x1 + (rho_1*x0 - rho_0*x1)*pow2_x2)/d;

      suma1 = suma1 + coef13*a*(pow(x2,3.f) - pow(x0,3.f)) + coef12*b*(pow2_x2 - pow2_x0) + c*(x2 - x0);
    }

 return suma1;
 }//function

