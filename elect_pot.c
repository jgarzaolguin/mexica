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

// Integral to evaluate the electrostatic potential for free atoms using GTOS basis 
// functions
//

double Elect_Pot_Free_GTO(int mu, int nu, int ang, int nt, double *matp, int *np, double *expo, double r)
{

double double_factorial(int);
double factorial_jul (int);
double exponential1 (int, double);
double constant_normalization_GTOS(int, double);
double numerical_incomplete_lower_gamma (double, double);
double incomplete_lower_gamma (double, float);
double incomplete_upper_gamma (double, float);
double suma_gamma (double, float);
double incomplete_analytical_lower_gamma (double, float);
double incomplete_analytical_upper_gamma (double, float);
double int_lower (double, int, int, int);
double int_upper (double, int, int, int);

int enes;
float k, p;
double arg2, inte1, inte2, Intesp, zetas, alpha;

zetas=expo[mu] + expo[nu]; 
enes=np[mu]+ np[nu];
alpha=zetas*pow(r,2);
        k=0.5*(enes + ang + 1);
        p=0.5*(enes - ang);
                inte1=pow(r,-ang-1)*0.5*pow((zetas),-k)*int_lower(alpha, np[mu], np[nu], ang);
                inte2=pow(r,ang)*0.5*pow(zetas,-p)*int_upper(alpha, np[mu], np[nu], ang);
                Intesp=constant_normalization_GTOS(np[mu], expo[mu])*constant_normalization_GTOS(np[nu], expo[nu])*(inte1+inte2);
//      printf("\nEl resultado de la integral espacial es de: %.12lf\n", Intesp);
return Intesp;
}

//Obtaining the value of the normalization constant for GTOS basis functions.

double doble_factorial (int ni)
{
 int I;
 double FAC;
 FAC=1.f;
for(I=1; I<=ni; I++)
{
        FAC*=(2.f*(double) I-1.f);
}
return FAC;
}

double exponential1(int ni, double expo)
{
        double RES;
        int e;
                e=ni+0.5;
                RES=pow(2*expo,e);
return RES;
}

double constant_normalization_GTOS (int ni, double expo)
{
        double RES, N;
                RES=pow(2,ni+1)*exponential1(ni, expo)/(doble_factorial(ni)*sqrt(M_PI));
                N=sqrt(RES);
return N;
}

double factorial_jul (int k)
{
 int I;
 double FAC;
 FAC=1.f;
 for (I=1; I<=k; I++)
   FAC*=(double)I;

 return FAC;
}

double numerical_incomplete_lower_gamma (double alpha, double z)
{
        double f1;
        double c0=0.189450610455068496285, c1=0.182603415044923588867, c2=0.169156519395002538189, c3=0.149595988816576732081, c4=0.124628971255533872052, c5=0.095158511682492784810, c6=0.062253523938647892863, c7=0.027152459411754094852, c8=c7, c9=c6, c10=c5, c11=c4, c12=c3, c13=c2, c14=c1, c15=c0;
        double x0=-0.095012509837637440185, x1=-0.281603550779258913230, x2=-0.458016777657227386342, x3=-0.617876244402643748447, x4=-0.755404408355003033895, x5=-0.865631202387831743880, x6=-0.944575023073232576078, x7=-0.989400934991649932596, x8=-x7, x9=-x6, x10=-x5, x11=-x4, x12=-x3, x13=-x2, x14=-x1, x15=-x0;
                f1=c0*pow(x0+1,z-1)*exp(-0.5*alpha*x0)+c1*pow(x1+1,z-1)*exp(-0.5*alpha*x1)+c2*pow(x2+1,z-1)*exp(-0.5*alpha*x2)+c3*pow(x3+1,z-1)*exp(-0.5*alpha*x3)+c4*pow(x4+1,z-1)*exp(-0.5*alpha*x4)+c5*pow(x5+1,z-1)*exp(-0.5*alpha*x5)+c6*pow(x6+1,z-1)*exp(-0.5*alpha*x6)+c7*pow(x7+1,z-1)*exp(-0.5*alpha*x7)+c8*pow(x8+1,z-1)*exp(-0.5*alpha*x8)+c9*pow(x9+1,z-1)*exp(-0.5*alpha*x9)+c10*pow(x10+1,z-1)*exp(-0.5*alpha*x10)+c11*pow(x11+1,z-1)*exp(-0.5*alpha*x11)+c12*pow(x12+1,z-1)*exp(-0.5*alpha*x12)+c13*pow(x13+1,z-1)*exp(-0.5*alpha*x13)+c14*pow(x14+1,z-1)*exp(-0.5*alpha*x14)+c15*pow(x15+1,z-1)*exp(-0.5*alpha*x15);
        return(f1);
}

double incomplete_lower_gamma (double alpha, float z)
{
        double lower_gamma;
                lower_gamma=pow(0.5*alpha,z)*exp(-0.5*alpha)*numerical_incomplete_lower_gamma (alpha ,z);
        return (lower_gamma);
}

double incomplete_upper_gamma (double alpha, float z)
{
        double upper_gamma, lower_gamma;
                lower_gamma=pow(0.5*alpha,z)*exp(-0.5*alpha)*numerical_incomplete_lower_gamma (alpha, z);
                upper_gamma=tgammaf(z)-lower_gamma;
        return (upper_gamma);
}

double suma_gamma (double alpha, float z)
{
        int I;
        double COS, RES;
for (I=0; I<z; I++)
{
        COS+=(pow(alpha,I)/factorial_jul(I));
}
        RES=exp(-alpha)*COS;
return RES;
}

double incomplete_analytical_lower_gamma (double alpha, float z)
{
        double lower_gamma;
                lower_gamma=tgamma(z)*(1-suma_gamma(alpha, z));
return lower_gamma;
}

double incomplete_analytical_upper_gamma (double alpha, float z)
{
        double upper_gamma;
                upper_gamma=tgamma(z)*suma_gamma(alpha, z);
return upper_gamma;
}

double int_lower (double alpha, int na, int nb, int l)
{
        int RES, x;
        float z;
        double int_1;
                x=na+nb+l;
                z=0.5*x+0.5;
                RES=pow(-1,x);
if(RES>0){
                int_1=incomplete_lower_gamma(alpha, z);
}else{
                int_1=incomplete_analytical_upper_gamma(alpha, z);
        }
return int_1;
}

double int_upper (double alpha, int na, int nb, int l)
{
        int RES, x;
        float z;
        double int_2;
                x=na+nb-l;
                z=0.5*x;
                RES=pow(-1,x);
        if(RES<0){
                int_2=incomplete_upper_gamma(alpha, z);
        }else{
                int_2=incomplete_analytical_upper_gamma(alpha, z);
                }
return int_2;
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
              double *NC_minus, double *NC_plus, char *basis)
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

  extern double Elect_Pot_Free_GTO(int mu, int nu, int ang, int nt, double *matp, int *np, double *expo, double r);

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
                  if (strcmp(basis,"STOs") == 0)
                    sumatot = sumatot + partial*Elect_Pot_Free_RHO(mu, nu, l, nt, matp, np, expo,
                                                             arreglo_factorial,
                                                             arreglo_inv_factorial, r);
                  else
		    sumatot= sumatot+partial*Elect_Pot_Free_GTO(mu, nu, l, nt, matp, np, expo, r);
                } else {
                if (strcmp(bound,"confined") == 0) {
                  sumatot = sumatot + partial*Elect_Pot_Impe_RHO(mu, nu, l, nt, matp, np,
                                                                 expo, arreglo_factorial,
                                                                 arreglo_inv_factorial, r, Rc);
                } else {
		if (strcmp(bound,"finite") == 0){
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
    }
  return -sumatot; // The sign is minus because we deal with the electron density
 }

int Evaluate_Elect_Pot(double z, int nt, double *matp, int *np, int *mang, int *ncm,
                       double *expo, char *bound, double *arreglo_factorial,
                       double *arreglo_inv_factorial, double *grid, int n_points, double Rc,
                       double *NC_minus, double *NC_plus, char *basis, double *pot_elect_grid)
{
  int i;
  double r, pot;
  extern double Elect_Pot_RHO(int nt, double *matp, int *np, int *mang, int *ncm,
                              double *expo, char *bound, double *arreglo_factorial,
                              double *arreglo_inv_factorial, double r, double Rc,
                              double *NC_minus, double *NC_plus, char *basis);

  for (i = 0; i < n_points; i++) {
    r = grid[i];
    pot = Elect_Pot_RHO(nt, matp, np, mang, ncm, expo, bound,
                        arreglo_factorial, arreglo_inv_factorial, r, Rc,
                        NC_minus, NC_plus, basis);
    pot_elect_grid[i] = pot;
  }

  return 0;
}

