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
// printf("I = %f \n", suma1);   // mike
 return suma1;
 }//function

double integral_two_points_interval(double x0, double x1, double f0, double f1) {
  double delta_x1x0, integral;
  
  delta_x1x0 = x1 - x0;
  integral = 0.5f*delta_x1x0*(f0 + f1);

  return(integral);
}

double integral_three_points_interval(double x0, double x1, double x2, double f0, double f1, double f2) {
  double total, num, denom, delta_x1x0, delta_x2x0, delta_x2x1;
  delta_x2x1 = x2 - x1;
  delta_x2x0 = x2 - x0;
  delta_x1x0 = x1 - x0;
  num = f1*delta_x2x0*delta_x2x0 - f0*delta_x2x1*(2.f*x0 - 3.f*x1 + x2) + f2*delta_x1x0*(x0 - 3.f*x1 + 2.f*x2);
  num = delta_x2x0*num;
  denom = 6.f*delta_x1x0*delta_x2x1;
  total = num/denom;
  return(total);
}

double integral_three_points(double *grid, double *function, int flag, int point_1, int point_2) {
  int i, bloques, indice, points, remain;
  double integral, x0, x1, x2, f0, f1, f2;

  if (flag == 1) points = point_1;
  else points = point_2;

  bloques = (points - 1)/2;
  integral = 0.f;
  for (i = 0; i < bloques; i++) {
    indice = 2*i;
    x0 = grid[indice];     f0 = function[indice];
    x1 = grid[indice + 1]; f1 = function[indice + 1];
    x2 = grid[indice + 2]; f2 = function[indice + 2];
    integral = integral + integral_three_points_interval(x0, x1, x2, f0, f1, f2);
  }

  remain = points - 1 - (indice + 2);

  if (remain == 1) {
    x0 = grid[indice + 2];  f0 = function[indice + 2];
    x1 = grid[indice + 3];  f1 = function[indice + 3];;
    integral = integral + integral_two_points_interval(x0, x1, f0, f1);
  }
//  printf("I3 = %f \n", integral);   // mike
  return(integral);
}
