#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//
//    Funciones para evaluar las integrales bielectronicas
//    Jorge Garza, Diciembre de 2006
//
//    Clebsch-Gordan function
//
      double clebsch_g(int j1, int j2, int j, int m1, int m2, int m)
      {
      int ent1, ent2, ent3, maximo1, maximo2, nup, ndown, mu, nexpo, ntest;
      double delta, rn1, rn2, rd1, rd2, term1, sum, adici, coef, total;

      extern int delta_kro_(int*, int*, double*);
      extern double factorial(int);

      ent1 = m;
      ent2 = m1 + m2;
      delta_kro_(&ent1, &ent2, &delta);
      ent1 = j + j2;
      ent2 = abs(j - j2);
      if (delta == (double)0 || j1 > ent1 || j1 < ent2)
        total = (double)0;
      else {
        rn1 = factorial(j - m)*factorial(j + m);
        rn1 = rn1*factorial(j1 + j2 - j)*factorial(j + j1 - j2);
        rn2 = factorial(j2 + j - j1);
        rd1 = factorial(j1 - m1)*factorial(j1 + m1);
        rd1 = rd1*factorial(j2 - m2)*factorial(j2 + m2);
        rd2 = factorial(j1 + j2 + j + 1);
        term1 = sqrt(rn1*rn2/(rd1*rd2));
        maximo1 = fmin(j + j2 + m1,j + m);
//      maximo1 = ( (j + j2 + m1) < (j + m) ? j + j2 + m1 : j + m);
          maximo2 = fmin(j + j2 + m1,j + j2 - j1);
//      maximo2 = ( (j + j2 + m1) > (j + j2 - j1) ? j + j2 + m1 : j + j2 - j1);
        nup = fmin(maximo1, maximo2);
//      nup = ( (maximo1, maximo2) ? maximo1 : maximo2);
        ndown = fmax(m1 - j1,m + j2 - j1);
//      ndown = ( (m1 - j1) > (m + j2 - j1) ? m1 - j1 : (m + j2 - j1));
        if (ndown < 0) ndown = 0;
        sum = (double)0;
        for (mu = ndown; mu <= nup; mu++) {
          rn1=factorial(j + j2 + m1 - mu)*factorial(j1 - m1 + mu);
          rd1=factorial(mu)*factorial(j + m - mu)*factorial(j + j2 - j1 - mu);
          rd2=factorial(mu + j1 - j2 - m);
          nexpo = mu + j2 + m2;
          if (nexpo != 0) {
            ntest = fmod(nexpo, 2);
            if (ntest == 0)
              coef = (double)1;
            else
              coef = -(double)1;
           }
          else
            coef = (double)1;
          adici = coef*rn1/(rd1*rd2);
          sum = sum + adici;
        }
        total = sqrt((double)(2*j + 1))*term1*sum;
      }
      return(total);
      }
//
      int delta_kro_(int* i, int* j, double* resultado)
      {
       int interno_i, interno_j;
       double total;
       interno_i = *i;
       interno_j = *j;
       total = (interno_i == interno_j ? (double)1 : (double)0);
       *resultado = total;

       return 1;
      }
//
//    Evaluacion de las integrales bielectronicas
//    Probablemente se puede evitar la multiplicacion de (2*l + 1)
//    en int3Y_1 y en int3Y_2, hay que intentar optimizarlo.
//
      double doselec(char   *using_gamma,
                     int     mu, 
                     int     nu, 
                     int     lam, 
                     int     sig, 
                     double  Rc, 
                     double  inv_Rc,
                     double  gamma_couple, 
                     char   *bound,
                     double *expo, 
                     int    *np, 
                     int    *ang, 
                     int    *ncm, 
                     double *zeta,
                     double *N_minus, 
                     double *N_plus,
                     double *arreglo_factorial, 
                     double *arreglo_inv_factorial)
      {
      int entero1, l, m, dif1, dif2, down, up, pot, lmu, lnu, llam, lsig;
      int mmu, mnu, mlam, msig, ndim, test;
      int sumam1, sumam2, compara1, compara2, compara3;
      double doble1, prod, sumatot, coef1, coef2, total, sumint, partrad, pi;

      extern double bielectronic_integral_confined(int, int, int, int, int, double, double, double*, int*, double*, double*);
      extern double bielectronic_integral_free(int, int, int, int, int, int*, double*, double*, double*);
      extern double int3Y_1(int, int, int, int, int ,int );
      extern double int3Y_2(int, int, int, int, int ,int );
      extern double bielectronic_integral_finite(char  *using_gamma,
                                                 int    ang, 
                                                 int    mu, 
                                                 int    nu, 
                                                 int    lam, 
                                                 int    sig, 
                                                 int    *np, 
                                                 int    *mang, 
                                                 double *zeta, 
                                                 double *alfa, 
                                                 double  Rc, 
                                                 double *arreglo_factorial, 
                                                 double *arreglo_inv_factorial,
                                                 double  gamma_couple,
                                                 double *N_minus, 
                                                 double *N_plus);

      doble1 = (double)1;
      pi = atan(doble1)*(double)4;
      lmu = ang[mu];
      lnu = ang[nu];
      llam = ang[lam];
      lsig = ang[sig];
      mmu = ncm[mu];
      mnu = ncm[nu];
      mlam = ncm[lam];
      msig = ncm[sig];
      dif1 = abs(lnu - lmu);
      dif2 = abs(lsig - llam);
      down = (dif1 < dif2 ? dif1 : dif2);
      dif1 = lnu + lmu;
      dif2 = lsig + llam;
      up = (dif1 > dif2 ? dif1 : dif2);
      sumatot = (double)0;
      for (l = down; l <= up; l++) {
        entero1 = lnu + lmu + l;
        test = fmod(entero1, 2);
        entero1 = llam + lsig + l;
        test = test*fmod(entero1,2);
        if (test == 0) {
          doble1 = (double)4;
          coef1 = doble1*pi/(double)(2*l + 1);
          sumint = (double)0;
          for (m = -l; m <= l; m++) {
            sumam1 = mmu + m;
            sumam2 = msig + m;
            if (mnu == sumam1 && mlam == sumam2) {
              coef2 = int3Y_1(lnu, lmu, l, mnu, mmu, m);
              coef2 = coef2*int3Y_2(llam, lsig, l, mlam, msig, m);
              if (coef2 != (double)0) {
//                compara1 = strcmp(bound,"confined");
                if (strcmp(bound,"confined") == 0) 
                  partrad = bielectronic_integral_confined(l, mu, nu, lam, sig, Rc, inv_Rc, expo, np, arreglo_factorial,
                                    arreglo_inv_factorial);
                else 
                if (strcmp(bound,"free") == 0) {
                  partrad = bielectronic_integral_free(l, mu, nu, lam, sig, np, expo,
                                    arreglo_factorial, arreglo_inv_factorial);
                }
                else
//CAMBIO EXPO POR ZETA
                 //partrad = bielectronic_integral_finite(l, mu, nu, lam, sig, np,
                 //                        ang, zeta, expo, Rc, 
                 //                        arreglo_factorial, arreglo_inv_factorial,
                 //                        gamma_couple,
                 //                        N_minus,  N_plus);
                  partrad = bielectronic_integral_finite(using_gamma,
                                                         l, 
                                                         mu, 
                                                         nu, 
                                                         lam, 
                                                         sig, 
                                                         np,
                                                         ang, 
                                                         expo, 
                                                         zeta, 
                                                         Rc, 
                                                         arreglo_factorial, 
                                                         arreglo_inv_factorial,
                                                         gamma_couple,
                                                         N_minus,  
                                                         N_plus);

              }
              else
                partrad = (double)0;
              coef2 = coef2*partrad;
            }
            else
              coef2 = (double)0;
            sumint = sumint + coef2;
          }
          total = coef1*sumint;
        }
        else
          total = (double)0;
        sumatot = sumatot + total;
      }
      return(sumatot);
      }
//
      double uno(int exp)
      {
       double total;
//     Tambien puede ser exp % 2
//       if (fmod(exp,2) == 0)
//     hay que verificarlo en un calculo

         total = (fmod(exp,2) == (double)0 ? (double)1 : -(double)1);
         return(total);
       }
//
//    La funcion int3Y evalua la integral de tres esfericos armonicos
//    Int[Y(l1,m1)Y(l2,m2)Y(l3,m3)]=
//    Sqrt[(2l1+1)(2l2+1)(2l3+1)/4Pi]*3j-Sym({l1,0},{l2,0},{l3,0})*
//    3j-Sym({l1,m1},{l2,m2},{l3,m3})
//    Para esto se usan los coeficientes de clebsch-Gordan.
//    La relacion que se esta usando es
//    coef({l1,m1},{l2,m2},{l3,m3})=
//    (-1)**(m3+l1-l2)Sqrt(2*l3+1)*3j-Sym({l1,m1},{l2,m2},{l3,m3})
//
//    Detalles sobre esta integral se pueden encontrar en
//    J. Phys. A:Math. Gen. vol. 31, pag. 7157 (1998)
//
      double int3Y(int l1, int l2, int l3, int m1, int m2, int m3)
      {
       int arg;
       double prod, pi;
       extern double clebsch_g(int, int, int, int, int, int);
       extern double uno(int);

       pi = ((double)4)*atan((double)1);
       prod = (double)((2*l1 + 1)*(2*l2 + 1));
       prod = prod/(double)(2*l3 + 1);
       prod = sqrt(prod/(((double)4)*pi));
       prod = prod*clebsch_g(l1, l2, l3, 0, 0, 0);
       prod = prod*clebsch_g(l1, l2, l3, m1, m2, -m3);
       arg = m3 - 2*l1 + 2*l2;
       prod = prod*uno(arg);
       return(prod);
      }
//
//    Se evalua la integral de tres esfericos armonicos estando
//    conjugados el segundo y el tercero
//    Int[Y(l1,m1)Y*(l2,m2)Y*(l3,m3)]=
//                 ((-1)**(m2+m3))int3Y((l1, l2, l3, m1, m2,m3)

      double int3Y_1(int l1, int l2, int l3, int m1, int m2, int m3)
      {
       int arg;
       double result;
       extern double int3Y(int, int, int, int, int, int);
       extern double uno(int);

       arg = m2 + m3;
       result = uno(arg)*int3Y(l1, l2, l3, m1, -m2, -m3);
       return(result);
       }
//
//    Se evalua la integral de tres esfericos armonicos estando
//    conjugado el primero
//    Int[Y*(l1,m1)Y(l2,m2)Y(l3,m3)]=
//                 ((-1)**m)int3Y((l1, l2, l3, m1, m2,m3)

      double int3Y_2(int l1, int l2, int l3, int m1, int m2, int m3)
      {
       int arg;
       double result;
       extern double int3Y(int, int, int, int, int, int);
       extern double uno(int);

       arg = m1;
       result = uno(arg)*int3Y(l1, l2, l3, -m1, m2, m3);
       return(result);
       }

