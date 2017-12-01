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
 double rho_1;
 double rho_2;
 double rho_3;
 double p1;
 double p2;
 double p3;
 double n;
 double l;
 double m;
 double powp1;
 double powp2;
 double powp3;
 double coef13;
 double coef43;
 int iglobal;

 coef13   = 1.f/3.f;
 coef43   = 4.f*coef13;

 suma1=0.f;

      for (i = 0; i < (points/2 - 1); i++) {
         iglobal = 2*(i + 1);

             p1 = grid[iglobal -2];
             p2 = grid[iglobal -1];
             p3 = grid[iglobal];

             rho_1 = grid_fun[iglobal -2];
             rho_2 = grid_fun[iglobal -1] ;
             rho_3 = grid_fun[iglobal];

      rho_1  = p1*p1*rho_1;
      rho_2  = p2*p2*rho_2;
      rho_3  = p3*p3*rho_3;

      powp1 = pow(p1,2.f);
      powp2 = pow(p2,2.f);
      powp3 = pow(p3,2.f);

        n=((p2 + p1)*(p3 + p2)/(p1*p2*(p3 + p2) - p2*p3*(p2 + p1)))*((powp2*rho_1 - powp1*rho_2)/(powp2 - powp1) - (powp3*rho_2 - powp2*rho_3)/(powp3 - powp2));

        l=(1.f/(p3*(p2 + p1) - p1*(p3 + p2)))*(p3*(powp2*rho_1 - powp1*rho_2)/(p2 - p1) - p1*(powp3*rho_2 - powp2*rho_3)/(p3 - p2));

        m=(rho_3 - n*p3 - l)/powp3;

        suma1 = suma1 + coef13*m*(pow(p3,3.f) - pow(p1,3.f)) + (1.f/2.f)*n*(powp3 - powp1) + l*(p3 - p1);
    }

 return suma1;
 }//function

