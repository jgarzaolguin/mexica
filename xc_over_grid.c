#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int get_varS(int points, double *rho, double *der_rho, double *varS) {
  int i;
  double one_3, pi, limit, rkf;
  one_3 = 1.f/3.f;
  pi = 4.f*atan(1.f);
  limit = 1e-20;
  for (i = 0; i < points; i++) {
    if (rho[i] <= limit) 
      varS[i] = 1e12;
    else {
      rkf = 3.f*pi*pi*rho[i];
      rkf = pow(rkf,one_3);
      varS[i] = fabs(der_rho[i])/(2.f*rkf*rho[i]);
    }
  }
  return 0;
}

double primder(int flag, int ini, double *grid, double *objetive) {
  double x0, x1, x2, x3, x4, y0, y1, y2, y3, y4, der;
  
  x0 = grid[ini]; 
  x1 = grid[ini + 1]; 
  x2 = grid[ini + 2]; 
  x3 = grid[ini + 3]; 
  x4 = grid[ini + 4]; 
  y0 = objetive[ini]; 
  y1 = objetive[ini + 1]; 
  y2 = objetive[ini + 2]; 
  y3 = objetive[ini + 3]; 
  y4 = objetive[ini + 4]; 

  switch(flag) {
  case 0:  der = y0/(x0 - x1) + y0/(x0 - x2) + y0/(x0 - x3) + y0/(x0 - x4) - 
           ((x0 - x2)*(x0 - x3)*(x0 - x4)*y1)/((x0 - x1)*(x1 - x2)*(x1 - x3)*(x1 - x4)) - 
           ((x0 - x1)*(x0 - x3)*(x0 - x4)*y2)/((x0 - x2)*(-x1 + x2)*(x2 - x3)*(x2 - x4)) - 
           ((x0 - x1)*(x0 - x2)*(x0 - x4)*y3)/((x0 - x3)*(-x1 + x3)*(-x2 + x3)*(x3 - x4)) - 
           ((x0 - x1)*(x0 - x2)*(x0 - x3)*y4)/((x0 - x4)*(-x1 + x4)*(-x2 + x4)*(-x3 + x4));
           break;
  case 1:  der = ((x1 - x2)*(x1 - x3)*(x1 - x4)*y0)/((x0 - x1)*(x0 - x2)*(x0 - x3)*(x0 - x4)) -
           y1/(x0 - x1) - ((-x0 + x1)*y1)/((x0 - x1)*(x1 - x2)) -
           ((-x0 + x1)*y1)/((x0 - x1)*(x1 - x3)) - ((-x0 + x1)*y1)/((x0 - x1)*(x1 - x4)) -
           ((-x0 + x1)*(x1 - x3)*(x1 - x4)*y2)/((x0 - x2)*(-x1 + x2)*(x2 - x3)*(x2 - x4)) -
           ((-x0 + x1)*(x1 - x2)*(x1 - x4)*y3)/((x0 - x3)*(-x1 + x3)*(-x2 + x3)*(x3 - x4)) -
           ((-x0 + x1)*(x1 - x2)*(x1 - x3)*y4)/((x0 - x4)*(-x1 + x4)*(-x2 + x4)*(-x3 + x4));
           break;
  case 2:  der = ((-x1 + x2)*(x2 - x3)*(x2 - x4)*y0)/((x0 - x1)*(x0 - x2)*(x0 - x3)*(x0 - x4)) - 
           ((-x0 + x2)*(x2 - x3)*(x2 - x4)*y1)/((x0 - x1)*(x1 - x2)*(x1 - x3)*(x1 - x4)) - 
           y2/(x0 - x2) - ((-x0 + x2)*y2)/((x0 - x2)*(-x1 + x2)) - 
           ((-x0 + x2)*y2)/((x0 - x2)*(x2 - x3)) - ((-x0 + x2)*y2)/((x0 - x2)*(x2 - x4)) - 
           ((-x0 + x2)*(-x1 + x2)*(x2 - x4)*y3)/((x0 - x3)*(-x1 + x3)*(-x2 + x3)*(x3 - x4)) - 
           ((-x0 + x2)*(-x1 + x2)*(x2 - x3)*y4)/((x0 - x4)*(-x1 + x4)*(-x2 + x4)*(-x3 + x4));
           break;
  case 3:  der = ((-x1 + x3)*(-x2 + x3)*(x3 - x4)*y0)/((x0 - x1)*(x0 - x2)*(x0 - x3)*(x0 - x4)) - 
           ((-x0 + x3)*(-x2 + x3)*(x3 - x4)*y1)/((x0 - x1)*(x1 - x2)*(x1 - x3)*(x1 - x4)) - 
           ((-x0 + x3)*(-x1 + x3)*(x3 - x4)*y2)/((x0 - x2)*(-x1 + x2)*(x2 - x3)*(x2 - x4)) - 
           y3/(x0 - x3) - ((-x0 + x3)*y3)/((x0 - x3)*(-x1 + x3)) - 
           ((-x0 + x3)*y3)/((x0 - x3)*(-x2 + x3)) - ((-x0 + x3)*y3)/((x0 - x3)*(x3 - x4)) - 
           ((-x0 + x3)*(-x1 + x3)*(-x2 + x3)*y4)/((x0 - x4)*(-x1 + x4)*(-x2 + x4)*(-x3 + x4));
           break;
  case 4 : der = ((-x1 + x4)*(-x2 + x4)*(-x3 + x4)*y0)/((x0 - x1)*(x0 - x2)*(x0 - x3)*(x0 - x4)) - 
           ((-x0 + x4)*(-x2 + x4)*(-x3 + x4)*y1)/((x0 - x1)*(x1 - x2)*(x1 - x3)*(x1 - x4)) - 
           ((-x0 + x4)*(-x1 + x4)*(-x3 + x4)*y2)/((x0 - x2)*(-x1 + x2)*(x2 - x3)*(x2 - x4)) - 
           ((-x0 + x4)*(-x1 + x4)*(-x2 + x4)*y3)/((x0 - x3)*(-x1 + x3)*(-x2 + x3)*(x3 - x4)) - 
           y4/(x0 - x4) - ((-x0 + x4)*y4)/((x0 - x4)*(-x1 + x4)) - 
           (-x0 + x4)*y4/((x0 - x4)*(-x2 + x4)) - (-x0 + x4)*y4/((x0 - x4)*(-x3 + x4));
           break;
  }

  return der;
}



int numerical_der(int points, double *grid, double *objetive, double *derivative) {
  int i, ini, flag;

  ini = 0;
  flag = 0;
  derivative[0] = primder(flag, ini, grid, objetive);
  flag = 1;
  derivative[1] = primder(flag, ini, grid, objetive);
  flag = 2;
  for (i = 0; i < (points - 4); i++) {
    derivative[i + 2] = primder(flag, i, grid, objetive); 
  } 
  ini = points - 5;
  flag = 3;
  derivative[points - 2] = primder(flag, ini, grid, objetive); 
  flag = 4;
  derivative[points - 1] = primder(flag, ini, grid, objetive); 

  return 0;
}

int xc_over_grid(int compara, char **save_dft, int flag_dft, double *weight_dft,  int points, double *grid, 
                 double *rho, double *der_rho, double *secder_rho) {
 int i, dft;
 double weight, pot_x, arg1, arg2, arg3, arg_derho, arg_grid, arg_secderho, arg_abs_grad_rho,
        arg_derabs_grad_rho;
 double varS[10000], dervarS[10000], abs_grad_rho[10000], derabs_grad_rho[10000];
 extern void pbegrid(double *rho, double *der_rho, double *arg_secderho,
                     double *abs_grad_rho, double *derabs_grad_rho,  double *varS,
                     double *der_varS, double *grid, double *pot_x);

 if (compara == 0) {
   for (dft = 1; dft < flag_dft; dft++) {
     weight = weight_dft[dft];
     if(strcmp(save_dft[dft], "pbe") == 0) {
       get_varS(points, rho, der_rho, varS);
       numerical_der(points, grid, varS, dervarS);
       for (i = 0; i < points; i++) abs_grad_rho[i] = fabs(der_rho[i]);
       numerical_der(points, grid, abs_grad_rho, derabs_grad_rho);
       for (i = 0; i < points; i++) {
         arg1 = rho[i];
         arg_derho = der_rho[i];
         arg_secderho = secder_rho[i];
         arg_abs_grad_rho = abs_grad_rho[i];
         arg_derabs_grad_rho = derabs_grad_rho[i];
         arg2 = varS[i];
         arg3 = dervarS[i];
         arg_grid = grid[i];
         pbegrid_(&arg1, &arg_derho, &arg_secderho, &arg_abs_grad_rho, &arg_derabs_grad_rho,
                  &arg2, &arg3, &arg_grid, &pot_x);
         printf("%f %f %f %f %f %f\n", grid[i], varS[i], dervarS[i], arg_secderho, derabs_grad_rho[i],
                 pot_x);
       }
     }
   } 
 } 
 else
   printf("No XC potential for open-shell atoms\n");

  return 0;
}
