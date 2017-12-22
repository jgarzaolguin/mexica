// Functions to evaluate the electrostatic potential on a grid.
// Jorge Garza, Dec/2017
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

// Getting Spherical Harmonics for \theta=0 and \phi=0
//
double Spherical_Harmonics_at_zero(int l, int m) {
  double total;
  if (l == 0) 
    total = 0.2820947917738781f;
  if (l == 1) {
    if (m == 0) total = 0.4886025119029199f;
    if (abs(m) == 1) total = 0.f;
  }
  if (l == 2) {
    if (m == 0) total = 0.6307831305050400f;
    if (abs(m) == 1) total = 0.f;
    if (abs(m) == 2) total = 0.f;
  }
  if (l == 3) {
    if (m == 0) total = 0.7463526651802308f;
    if (abs(m) == 1) total = 0.f;
    if (abs(m) == 2) total = 0.f;
    if (abs(m) == 3) total = 0.f;
  }
  if (l == 4) {
    if (m == 0) total = 0.8462843753216344f;
    if (abs(m) == 1) total = 0.f;
    if (abs(m) == 2) total = 0.f;
    if (abs(m) == 3) total = 0.f;
    if (abs(m) == 4) total = 0.f;
  }
  if (l == 5) {
    if (m == 0) total = 0.9356025796273888f;
    if (abs(m) == 1) total = 0.f;
    if (abs(m) == 2) total = 0.f;
    if (abs(m) == 3) total = 0.f;
    if (abs(m) == 4) total = 0.f;
    if (abs(m) == 5) total = 0.f;
  }
  if (l == 6) {
    if (m == 0) total = 1.017107236282055f;
    if (abs(m) == 1) total = 0.f;
    if (abs(m) == 2) total = 0.f;
    if (abs(m) == 3) total = 0.f;
    if (abs(m) == 4) total = 0.f;
    if (abs(m) == 5) total = 0.f;
    if (abs(m) == 6) total = 0.f;
  }
  if (l == 7) {
    if (m == 0) total = 1.092548430592079f;
    if (abs(m) == 1) total = 0.f;
    if (abs(m) == 2) total = 0.f;
    if (abs(m) == 3) total = 0.f;
    if (abs(m) == 4) total = 0.f;
    if (abs(m) == 5) total = 0.f;
    if (abs(m) == 6) total = 0.f;
    if (abs(m) == 7) total = 0.f;
  }
  if (l == 8) {
    if (m == 0) total = 1.163106622920320f;
    if (abs(m) == 1) total = 0.f;
    if (abs(m) == 2) total = 0.f;
    if (abs(m) == 3) total = 0.f;
    if (abs(m) == 4) total = 0.f;
    if (abs(m) == 5) total = 0.f;
    if (abs(m) == 6) total = 0.f;
    if (abs(m) == 7) total = 0.f;
    if (abs(m) == 8) total = 0.f;
  }
  return total;
}

// Integral to evaluate the electrostatic potential
// for free atoms.
//
double Elect_Pot_Free_RHO(int mu, int nu, int ang, int nt, double *matp, int *np,
                      double *expo, double *arreglo_factorial,
                      double *arreglo_inv_factorial, double r)
{
  int arg1;
  double arg2, result1, result2, total, int_ep_1, int_ep_2, int_ep_3, suma;
  extern int constant_normalization_free(int*, double*, int*, double*);
  extern double int1l(int arg1, double arg2, double *arreglo_factorial);
  extern double intc(int a, double b, double r, double *arreglo_factorial,
             double *arreglo_inv_factorial);

  constant_normalization_free(&mu, expo, np, &result1);
  constant_normalization_free(&nu, expo, np, &result2);
// Debug      printf("mu, nu = %d, %d; const = %f, %f\n", mu, nu, result1, result2);
// Debug      printf("n_mu, n_nu = %d, %d; e_mu, e_nu = %f, %f\n", np[mu], np[nu], expo[mu], expo[nu]);
  total = result1*result2;

  arg2 = expo[mu] + expo[nu];

  arg1 = np[mu] + np[nu] - ang - 1;
  // \int_0^{\infty} dx x^{n_{\mu} + n_{\nu} - l - 1} e^{-(\xi_{\mu} + \xi_{\nu})*x}
  int_ep_2 = int1l(arg1, arg2, arreglo_factorial); 

  if (r == 0.f) {
    if (ang == 0) total = total*int_ep_2;
    else total = 0.f;
  } else {
      arg1 = np[mu] + np[nu] + ang;
      // \int_0^r dx x^{n_{\mu} + n_{\nu} + l} e^{-(\xi_{\mu} + \xi_{\nu})*x}
      int_ep_1 = intc(arg1, arg2, r, arreglo_factorial, arreglo_inv_factorial); 
      arg1 = np[mu] + np[nu] - ang - 1;
      // \int_0^r dx x^{n_{\mu} + n_{\nu} - l - 1} e^{-(\xi_{\mu} + \xi_{\nu})*x}
      int_ep_3 = intc(arg1, arg2, r, arreglo_factorial, arreglo_inv_factorial); 
      if (ang == 0)
        total = total*(int_ep_1/r + int_ep_2 - int_ep_3);
      else
        total = total*(int_ep_1/pow(r,(double) ang + 1) + 
                     pow(r,(double) ang)*(int_ep_2 - int_ep_3));
  }

  return total;
 }

// Integral to evaluate the electrostatic potential
// for atoms confined by impenetrable walls.
//
double Elect_Pot_Impe_RHO(int mu, int nu, int ang, int nt, double *matp, int *np,
                      double *expo, double *arreglo_factorial,
                      double *arreglo_inv_factorial, double r, double Rc)
{
  int arg1;
  double arg2, result1, result2, total, int_ep_1_1, int_ep_2_1, int_ep_3_1,
         int_ep_1_2, int_ep_2_2, int_ep_3_2, int_ep_1_3, int_ep_2_3, int_ep_3_3,
         int_ep_1, int_ep_2, int_ep_3, suma,
         inv_Rc, inv_Rc_cuad;
  extern int constant_normalization_free(int*, double*, int*, double*);
  extern double intc(int a, double b, double r, double *arreglo_factorial,
             double *arreglo_inv_factorial);

  constant_normalization_free(&mu, expo, np, &result1);
  constant_normalization_free(&nu, expo, np, &result2);
// Debug      printf("mu, nu = %d, %d; const = %f, %f\n", mu, nu, result1, result2);
// Debug      printf("n_mu, n_nu = %d, %d; e_mu, e_nu = %f, %f\n", np[mu], np[nu], expo[mu], expo[nu]);
  total = result1*result2;

  arg2 = expo[mu] + expo[nu];

  inv_Rc = 1.f/Rc;
  inv_Rc_cuad = inv_Rc/Rc;
  arg1 = np[mu] + np[nu] - ang - 1;
  // \int_0^{R_c} dx x^{n_{\mu} + n_{\nu} - l - 1} e^{-(\xi_{\mu} + \xi_{\nu})*x}
  int_ep_1_2 = intc(arg1, arg2, Rc, arreglo_factorial, arreglo_inv_factorial);
  arg1 = np[mu] + np[nu] - ang;
  // \int_0^{Rc} dx x^{n_{\mu} + n_{\nu} - l} e^{-(\xi_{\mu} + \xi_{\nu})*x}
  int_ep_2_2 = intc(arg1, arg2, Rc, arreglo_factorial, arreglo_inv_factorial);
  arg1 = np[mu] + np[nu] - ang + 1;
  // \int_0^{Rc} dx x^{n_{\mu} + n_{\nu} - l + 1} e^{-(\xi_{\mu} + \xi_{\nu})*x}
  int_ep_3_2 = intc(arg1, arg2, Rc, arreglo_factorial, arreglo_inv_factorial);
  int_ep_2 = int_ep_1_2 - 2.f*inv_Rc*int_ep_2_2 + inv_Rc_cuad*int_ep_3_2;

  if (r == 0.f) {
    if (ang == 0) total = total*int_ep_2;
    else total = 0.f;
  } else {
      arg1 = np[mu] + np[nu] + ang;
      // \int_0^r dx x^{n_{\mu} + n_{\nu} + l} e^{-(\xi_{\mu} + \xi_{\nu})*x}
      int_ep_1_1 = intc(arg1, arg2, r, arreglo_factorial, arreglo_inv_factorial);
      arg1 = np[mu] + np[nu] + ang + 1;
      // \int_0^r dx x^{n_{\mu} + n_{\nu} + l + 1} e^{-(\xi_{\mu} + \xi_{\nu})*x}
      int_ep_2_1 = intc(arg1, arg2, r, arreglo_factorial, arreglo_inv_factorial);
      arg1 = np[mu] + np[nu] + ang + 2;
      // \int_0^r dx x^{n_{\mu} + n_{\nu} + l + 2} e^{-(\xi_{\mu} + \xi_{\nu})*x}
      int_ep_3_1 = intc(arg1, arg2, r, arreglo_factorial, arreglo_inv_factorial);
      int_ep_1 = int_ep_1_1 - 2.f*inv_Rc*int_ep_2_1 + inv_Rc_cuad*int_ep_3_1;
      arg1 = np[mu] + np[nu] - ang - 1;
      // \int_0^r dx x^{n_{\mu} + n_{\nu} - l - 1} e^{-(\xi_{\mu} + \xi_{\nu})*x}
      int_ep_1_3 = intc(arg1, arg2, r, arreglo_factorial, arreglo_inv_factorial);
      arg1 = np[mu] + np[nu] - ang;
      // \int_0^r dx x^{n_{\mu} + n_{\nu} - l} e^{-(\xi_{\mu} + \xi_{\nu})*x}
      int_ep_2_3 = intc(arg1, arg2, r, arreglo_factorial, arreglo_inv_factorial);
      arg1 = np[mu] + np[nu] - ang + 1;
      // \int_0^r dx x^{n_{\mu} + n_{\nu} - l + 1} e^{-(\xi_{\mu} + \xi_{\nu})*x}
      int_ep_3_3 = intc(arg1, arg2, r, arreglo_factorial, arreglo_inv_factorial);
      int_ep_3 = int_ep_1_3 - 2.f*inv_Rc*int_ep_2_3 + inv_Rc_cuad*int_ep_3_3;
      if (ang == 0)
        total = total*(int_ep_1/r + int_ep_2 - int_ep_3);
      else
        total = total*(int_ep_1/pow(r,(double) ang + 1) +
                     pow(r,(double) ang)*(int_ep_2 - int_ep_3));
  }

  return total;
 }

// Integral to evaluate the electrostatic potential
// for atoms confined by impenetrable walls.
//
double Elect_Pot_Pen_RHO(int mu, int nu, int ang, int nt, double *matp, int *np,
                         int *mang, double *expo, double *arreglo_factorial,
                         double *arreglo_inv_factorial, double r, double Rc,
                         double *NC_minus, double *NC_plus)
{
  int enes, eles;
  double alphas, zetas, product_N_minus, product_N_plus, arg1, arg2, total,
         int_1, int_2, int_complete, part1, part2, part3, part4;

  extern double intc(int, double, double, double*, double*);
  extern double upper_incomplete_gamma(double, int, double);

  zetas = expo[mu] + expo[nu];
  enes = np[mu] + np[nu];
  eles = mang[mu] + mang[nu];
  product_N_minus = NC_minus[mu]*NC_minus[nu];
  product_N_plus  = NC_plus[mu]*NC_plus[nu];
  alphas = zetas - (double) (enes + eles)/Rc;
  arg1 = enes - ang - 1;
  arg2 = eles + ang + 1;
  int_complete = product_N_minus*intc(arg1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)
                 + product_N_plus*upper_incomplete_gamma(Rc, arg2, alphas);
  if (r <= Rc) {
    if (r == 0.f) {
      if (ang == 0)
        total = int_complete;
      else
        total = 0.f;
    } else {
      arg1 = enes + ang;
      int_1 = intc(arg1, zetas, r, arreglo_factorial, arreglo_inv_factorial);
      arg1 = enes - ang - 1;
      int_2 = intc(arg1, zetas, r, arreglo_factorial, arreglo_inv_factorial);
      if (ang == 0) 
        total = product_N_minus*(int_1/r - int_2) + int_complete;
      else 
        total = product_N_minus*(int_1/pow(r,(double) ang + 1) - pow(r,(double) ang)*int_2) +
                pow(r, (double) ang)*int_complete;
    } // Finish condition (r > 0 && r <= R_c)
  } // Finish condition r <= R_c
  else {
// This is for the region r > R_c
    arg1 = enes + ang;
    part1 = product_N_minus*intc(arg1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
    arg2 = eles - ang;
    int_1 = upper_incomplete_gamma(Rc, arg2, alphas);
    int_2 = upper_incomplete_gamma(r, arg2, alphas);
    part2 = product_N_plus*(int_1 - int_2);
    arg1 = enes - ang - 1;
    part3 = product_N_minus*intc(arg1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
    arg2 = eles + ang + 1;
    int_1 = upper_incomplete_gamma(Rc, arg2, alphas);
    int_2 = upper_incomplete_gamma(r, arg2, alphas);
    part4 = product_N_plus*(int_1 - int_2);
    if (ang == 0) {
      total = (part1 + part2)/r;
      total = total + int_complete;
      total = total - (part3 + part4);
    } else {
      total = (part1 + part2)/pow(r, (double) ang + 1);
      total = total + pow(r,(double) ang)*int_complete;
      total = total - pow(r,(double) ang)*(part3 + part4);
    }
  }
  return total;
 }

double Elect_Pot_RHO(int nt, double *matp, int *np, int *ang, int *ncm,
              double *expo, char *bound, double *arreglo_factorial,
              double *arreglo_inv_factorial, double r, double Rc,
              double *NC_minus, double *NC_plus)
{
  int l, m, mu, nu, lmu, lnu, mmu, mnu, down, up, entero1, test,
      sumam1, sumam2;
  double pi, doble1, coef1, coef2, sumint, sumatot, partial;

  extern double int3Y_1(int, int, int, int, int ,int );
  extern double int3Y_2(int, int, int, int, int ,int );
  extern double Elect_Pot_Free_RHO(int mu, int nu, int ang, int nt, double *matp, int *np,
                      double *expo, double *arreglo_factorial,
                      double *arreglo_inv_factorial, double r);
  extern double Elect_Pot_Impe_RHO(int mu, int nu, int ang, int nt, double *matp, int *np,
                      double *expo, double *arreglo_factorial,
                      double *arreglo_inv_factorial, double r, double Rc);
  extern double Elect_Pot_Pen_RHO(int mu, int nu, int l, int nt, double *matp, int *np,
                         int *ang, double *expo, double *arreglo_factorial,
                         double *arreglo_inv_factorial, double r, double Rc,
                         double *NC_minus, double *NC_plus);

  doble1 = (double)1;
  pi = atan(doble1)*(double)4;
  
  sumatot = (double)0;
  for (mu = 0; mu < nt; mu++)
    for (nu = 0; nu < nt; nu++) {
      lmu = ang[mu];
      lnu = ang[nu];
      mmu = ncm[mu];
      mnu = ncm[nu];
      down = abs(lnu - lmu);
      up = lnu + lmu;
      for (l = down; l <= up; l++) {
        entero1 = lnu + lmu + l;
        test = fmod(entero1, 2);
        if (test == 0) {
          doble1 = (double)4;
          coef1 = doble1*pi/((double)(2*l + 1));
          sumint = (double)0;
          for (m = -l; m <= l; m++) {
            sumam1 = mmu + m;
            if (mnu == sumam1) {
              coef2 = int3Y_1(lnu, lmu, l, mnu, mmu, m);
              if (coef2 != (double)0) {
                partial = matp[mu*nt + nu]*coef1*coef2*Spherical_Harmonics_at_zero(l,m);
                if (strcmp(bound,"free") == 0) {
                  sumatot = sumatot + partial*Elect_Pot_Free_RHO(mu, nu, l, nt, matp, np, expo,
                                                             arreglo_factorial,
                                                             arreglo_inv_factorial, r);
                } else {
                if (strcmp(bound,"confined") == 0) {
                  sumatot = sumatot + partial*Elect_Pot_Impe_RHO(mu, nu, l, nt, matp, np,
                                                                 expo, arreglo_factorial,
                                                                 arreglo_inv_factorial, r, Rc);
                } else {
                  sumatot = sumatot + partial*Elect_Pot_Pen_RHO(mu, nu, l, nt, matp, np, ang, expo,
                                                        arreglo_factorial,
                                                        arreglo_inv_factorial, r, Rc,
                                                        NC_minus, NC_plus);
                }
                }
              }
            }
          }
        }
      }
    }
  return -sumatot; // The sign is minus because we deal with the electron density
 }

int Evaluate_Elect_Pot(double z, int nt, double *matp, int *np, int *mang, int *ncm,
                       double *expo, char *bound, double *arreglo_factorial,
                       double *arreglo_inv_factorial, double *grid, int n_points, double Rc,
                       double *NC_minus, double *NC_plus)
{
  int i;
  double r, pot;
  extern double Elect_Pot_RHO(int nt, double *matp, int *np, int *mang, int *ncm,
                              double *expo, char *bound, double *arreglo_factorial,
                              double *arreglo_inv_factorial, double r, double Rc,
                              double *NC_minus, double *NC_plus);

  FILE *target;

  target = fopen("Elect_pot.dat","w");
  for (i = 0; i < n_points; i++) {
    r = grid[i];
    pot = Elect_Pot_RHO(nt, matp, np, mang, ncm, expo, bound,
                        arreglo_factorial, arreglo_inv_factorial, r, Rc,
                        NC_minus, NC_plus);
    fprintf(target, "%20.8f  %20.8f\n", grid[i], pot/2.f);
  }

  fclose(target);

  return 0;
}

