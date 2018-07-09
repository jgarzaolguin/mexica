#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//IMPENETRABLE STO(1-r/RC)


 double msto(int     mu, 
             double  r, 
             double  Rc, 
             double *expo, 
             int    *np, 
             double *arreglo_factorial,
             double *arreglo_inv_factorial)
 {
  double total, constant;

  extern double constc(int     mu, 
                       double  Rc, 
                       double *expo, 
                       int    *np, 
                       double *arreglo_factorial,
                       double *arreglo_inv_factorial);
  constant = constc(mu, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);


  if (r == 0) {
    if (np[mu] == 1)
     total = constant;
    else
     total = 0.f;
  } else { 
      total = constant*(1.f - r/Rc)*exp(-expo[mu]*r); // total = N[mu]*exp(-expo[mu]*r)*(1-r/Rc)
       if (np[mu] == 1)
         total = total;
       else
        total = total*pow(r,(double) (np[mu] - 1));   //Aqui total = N[mu]*r^(np[mu] - 1)* exp(-expo[mu]*r)
  }

  return (total);
}
///////////////////////////////////////////////////
//Primera derivada del msto
 double der_msto_r(int mu, double r, double Rc, double* expo, int* np, double *arreglo_factorial,
                   double *arreglo_inv_factorial)
 {
  double total, factor, constant;
  extern double constc(int mu, double Rc, double* expo, int* np, double *arreglo_factorial,
                       double *arreglo_inv_factorial);
  constant = constc(mu, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);

  if (r == 0.f) {
    if (np[mu] == 1)
      total = -constant*(1.f/Rc + expo[mu]);
    else
    if (np[mu] == 2)
      total = constant;
    else
      total = 0.f;
  } else {

    total = constant*exp(-expo[mu]*r);
    factor = (1.f - r/Rc)*expo[mu] + 1.f/Rc;

    if (np[mu] == 1) 
      total = -total*factor;
   
    else
    if (np[mu] == 2) {
      total = total*(-2.f*r + Rc + r*(r - Rc)*expo[mu]);
      total = total/Rc;
    }
    else {
      total = total*pow(r,(double) (np[mu] - 2));
      total = total*(Rc + (r - Rc)*((double) np[mu] - r*expo[mu]));
      total = total/Rc;
    }
  } 
  return (total); 
  }
//Segunda derivada de la función msto
 double sec_der_msto_r(int mu, double r, double Rc, double* expo, int* np, double *arreglo_factorial,
                       double *arreglo_inv_factorial)
 {
  double total, factor, constant;

  extern double constc(int mu, double Rc, double* expo, int* np, double *arreglo_factorial,
                       double *arreglo_inv_factorial);
  constant = constc(mu, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
  if (r == 0.f) {
    if (np[mu] == 1) {
      total = constant*expo[mu]*(2.f + Rc*expo[mu]);
      total = total/Rc; 
    } else
        if (np[mu] == 2) {
          total = -2.f*constant*(1.f + Rc*expo[mu]);
          total = total/Rc;
        }
         else
          if (np[mu] == 3)
            total = 2.f*constant;
             else
               total = 0.f;
  } else {
       total = exp(-r*expo[mu])*constant;
       if (np[mu] == 1) {
         total = total*expo[mu]*(2.f - r*expo[mu] + Rc*expo[mu]);
         total = total/Rc;
       } else
          if (np[mu] == 2) {
            total = total*(-pow(r*expo[mu],2.f) - 2.f*(1.f + Rc*expo[mu]) + r*expo[mu]*(4.f + Rc*expo[mu]));
            total = total/Rc;
          } else 
           if (np[mu] == 3) {
             total = total*(2.f*Rc - pow(r,3.f)*pow(expo[mu],2.f) + pow(r,2.f)*expo[mu]*(6.f + Rc*expo[mu]) - 2.f*r*(3.f + 2.f*Rc*expo[mu]));
             total = total/Rc;
           } else  {
              total  = total*pow(r,(double) (np[mu] - 3.f));
              factor = (double) pow((double) np[mu],2.f)*(r - Rc)
                      +(double) pow(r,3.f)*pow(expo[mu],2.f)
                      -(double) np[mu]*(r - 3.f*Rc + 2.f*pow(r,2.f)*expo[mu] - 2.f*r*Rc*expo[mu]) 
                      -(double) Rc*(2.f + 2.f*r+expo[mu] + pow(r*expo[mu],2.f));
              total = total*factor;
             }
    }

  return (total);
  }
//Orbital confinado
 double confined_orbital(int nt, int orbital, double r, double Rc, double* expo,   
                         int* np, double* vectors, double *arreglo_factorial,      
                         double *arreglo_inv_factorial)                            
 {
  int i;
  double suma;
  extern double msto(int mu, double r, double Rc, double* expo, int* n,
                     double *arreglo_factorial, double *arreglo_inv_factorialp);

  suma = 0.f;
  for (i = 0; i < nt; i++) 
    suma = suma + vectors[i*nt + orbital]*msto(i, r, Rc, expo, np, arreglo_factorial,
                                               arreglo_inv_factorial);

  return (suma);
  }

//Derivada del orbital confinado
 double der_confined_orbital(int nt, int orbital, double r, double Rc, double* expo,
                             int* np, double* vectors, double *arreglo_factorial,
                             double *arreglo_inv_factorial)
 {
  int i;
  double suma;
  extern double der_msto_r(int mu, double r, double Rc, double* expo, int* n,
                           double *arreglo_factorial, double *arreglo_inv_factorial);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*der_msto_r(i, r, Rc, expo, np, arreglo_factorial,
                                                     arreglo_inv_factorial);
  }

  return (suma);
  }
// Segunda derivada del orbital
 double sec_der_confined_orbital(int nt, int orbital, double r, double Rc, double* expo,
                                 int* np, double* vectors, double *arreglo_factorial,
                                 double *arreglo_inv_factorial)
 {
  int i;
  double suma;
  extern double sec_der_msto_r(int mu, double r, double Rc, double* expo, int* n,
                               double *arreglo_factorial, double *arreglo_inv_factorial);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*sec_der_msto_r(i, r, Rc, expo, np, arreglo_factorial,
                                                         arreglo_inv_factorial);
  }

  return (suma);
  }
//Densidad radial
 double rho_radial_confined(int nt, int elec, double r, double Rc, double* expo,
                            int* np, double* vectors, char* tipo, double *arreglo_factorial,
                            double *arreglo_inv_factorial)
 {
  int orbital;
  double sum, partial, occ, pi;
  extern double confined_orbital(int nt, int orbital, double r, double Rc, double* expo,
                                 int* np, double* vectors, double *arreglo_factorial,
                                 double *arreglo_inv_factorial);

  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {
    partial = confined_orbital(nt, orbital, r, Rc, expo, np, vectors, arreglo_factorial,
                               arreglo_inv_factorial);
    sum = sum + occ*partial*partial;
  }

  return(sum/(4.f*pi));
 }
//Derivada de la densidad radial
 double der_rho_radial_confined(int nt, int elec, double r, double Rc, double* expo,
                                int* np, double* vectors, char* tipo, double *arreglo_factorial,
                                double *arreglo_inv_factorial)
 {
  int orbital;
  double sum, partial, occ, pi;
  extern double confined_orbital(int nt, int orbital, double r, double Rc, double* expo,
                                 int* np, double* vectors, double *arreglo_factorial,
                                 double *arreglo_inv_factorial);

  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {
    partial = 2.f*confined_orbital(nt, orbital, r, Rc, expo, np, vectors, arreglo_factorial,
                                   arreglo_inv_factorial);
    partial = partial*der_confined_orbital(nt, orbital, r, Rc, expo, np, vectors, arreglo_factorial,
                                           arreglo_inv_factorial);
    sum = sum + occ*partial;
  }

  return(sum/(4.f*pi));
 }

//Segunda derivada de la densidad radial
 double sec_der_rho_radial_confined(int nt, int elec, double r, double Rc, double* expo,
                                    int* np, double* vectors, char* tipo, double *arreglo_factorial,
                                    double *arreglo_inv_factorial)
 {
  int orbital;
  double sum, partial1, partial2, partial, occ, pi;
  extern double confined_orbital(int nt, int orbital, double r, double Rc, double* expo,
                                 int* np, double* vectors, double *arreglo_factorial,
                                 double *arreglo_inv_factorial);

  extern double sec_der_confined_orbital(int nt, int orbital, double r, double Rc, double* expo,
                                         int* np, double* vectors, double *arreglo_factorial,
                                         double *arreglo_inv_factorial);

  extern double der_confined_orbital(int nt, int orbital, double r, double Rc, double* expo,
                                     int* np, double* vectors, double *arreglo_factorial,
                                     double *arreglo_inv_factorial);
  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {
    partial1 = der_confined_orbital(nt, orbital, r, Rc, expo, np, vectors, arreglo_factorial,
                                    arreglo_inv_factorial);
  
    partial1 = pow(partial1, 2.f);
   
    partial2 = sec_der_confined_orbital(nt, orbital, r, Rc, expo, np, vectors, arreglo_factorial,
                                        arreglo_inv_factorial);

    partial2 = partial2*confined_orbital(nt, orbital, r, Rc, expo, np, vectors, arreglo_factorial,
                                         arreglo_inv_factorial);

    partial = 2.f*(partial1 + partial2);

    sum = sum + occ*partial;
  }

  return(sum/(4.f*pi));
 }
 /* Here I'm going to start with the GTOs, this case is only for the impenetrable wall */
 /* The R_{GTO}(r) = N^{(imp)}*r^{n - 1}*exp^{-\alpha r^{2}}*(1 - r/r0) */ 
 double gto_imp(int mu, double r, double r0, double* expo, int* np) {
    extern double cte_norm_gto_imp(int *, int *, double *, double , double *);
    double result, n_gto_imp, factor, num;

    cte_norm_gto_imp(&mu, np, expo, r0, &n_gto_imp);

    if(r == 0) {
       if(np[mu] == 1){
          result = n_gto_imp;
       }
       else{
          result = 0.f;
       }
    }
    else {
       factor = n_gto_imp*exp(-expo[mu]*r*r)*(1.f - r/r0);
       if(np[mu] == 1){
          result = factor;
       }
       else{
          num = (double)np[mu];
          result = pow(r,num - 1.f)*factor;
       }
    }
    return (result);
 }

 double der_gto_imp_r(int mu, double r, double r0, double* expo, int* np) {
    extern double cte_norm_gto_imp(int *, int *, double *, double , double *);
    double result, n_gto_imp, factor, num;

    cte_norm_gto_imp(&mu, np, expo, r0, &n_gto_imp);

    if(r == 0) {
       if(np[mu] == 1){
          result = -n_gto_imp/r0;
       }
       else{
          if(np[mu] == 2){
             result = n_gto_imp;
          }
          else{
             result = 0.f;
          }
       }
    }
    else{
       num = (double)np[mu];
       factor = n_gto_imp*exp(-expo[mu]*r*r);
       if(np[mu] == 1){
          result =  (-2.f*expo[mu]*r*(1.f - r/r0) - 1.f/r0)*factor;
       }
       else{
          if(np[mu] == 2){
             result = ((1.f - 2.f*expo[mu]*r*r)*(1.f - r/r0) - r/r0)*factor;
          }
          else{
             result = ((num - 1.f)*pow(r,num - 2.f) - 2.f*expo[mu]*pow(r,num))*(1.f - r/r0) - pow(r,num - 1.f)/r0;
             result = result*factor;
          }
       }
    }
  return (result);
  }

  double sec_der_gto_imp_r(int mu, double r, double r0, double* expo, int* np) {
    extern double cte_norm_gto_imp(int *, int *, double *, double , double *);
    double result, n_gto_imp, factor, num;

    cte_norm_gto_imp(&mu, np, expo, r0, &n_gto_imp);

    if(r == 0) {
       if(np[mu] == 1){
          result = -2.f*n_gto_imp*expo[mu];
       }
       else{
          if(np[mu] == 2){
             result = 2.f*n_gto_imp/r0;
          }
          else{
             if(np[mu] == 3){
                result = 2.f*n_gto_imp;   
             } 
             else{
                result = 0.f;
             } 
          }
       }
    }
    else{
       factor = n_gto_imp*exp(-expo[mu]*r*r);
       num = (double)np[mu];
       if(np[mu] == 1){
          result = 2.f*factor*expo[mu]*((2.f*expo[mu]*r*r - 1.f)*(1.f - r/r0) + 2.f*r/r0);
       }
       else{
          if(np[mu] == 2){
             result = 2.f*factor*(2.f*pow(expo[mu],2.f)*pow(r,3.f)*(1.f - r/r0) - (1.f - 2.f*expo[mu]*pow(r,2.f))/r0); 
          }
          else{
             if(np[mu] == 3){
                result = 2.f*factor*((1.f - 5.f*expo[mu]*pow(r,2.f) + 2.f*pow(expo[mu],2.f)*pow(r,4.f))*(1.f - r/r0) - (2.f*r/r0)*(1.f - expo[mu]*pow(r,2.f))); 
             }
             else{
                result = factor*pow(r,num - 3.f);
                result = result*(
                         ((num - 1.f)*(num - 2.f) + 2.f*expo[mu]*(1.f - 2.f*num)*pow(r,2.f) + 4.f*pow(expo[mu],2.f)*pow(r,4.f))*(1.f - r/r0) - 
                         (2.f*r/r0)*(num - 1.f - 2.f*expo[mu]*pow(r,2.f))  
                         );
             }
          }
       }
    }
  return (result);
  }

 double orbital_gto_imp(int nt, int orbital, double r, double r0, double* expo, int* np, double* vectors) {
    int i;
    double suma;
    extern double gto_imp(int , double , double , double* , int* );

    suma = 0.f;
    for(i = 0; i < nt; i++) {
       suma = suma + vectors[i*nt + orbital]*gto_imp(i, r, r0, expo, np);
    }

    return (suma);
 }

 double der_orbital_gto_imp(int nt, int orbital, double r, double r0, double* expo, int* np, double* vectors) {
    int i;
    double suma;
    extern double der_gto_imp_r(int , double , double , double* , int* ); 

    suma = 0.f;
    for(i = 0; i < nt; i++) {
       suma = suma + vectors[i*nt + orbital]*der_gto_imp_r(i, r, r0, expo, np);
    }

    return (suma);
 }

 double sec_der_orbital_gto_imp(int nt, int orbital, double r, double r0, double* expo, int* np, double* vectors) {
    int i;
    double suma;
    extern double sec_der_gto_imp_r(int , double , double , double* , int* );

    suma = 0.f;
    for(i = 0; i < nt; i++) {
       suma = suma + vectors[i*nt + orbital]*sec_der_gto_imp_r(i, r, r0, expo, np);
    }

    return (suma);
 }

 double rho_radial_gto_imp(int nt, int elec, double r, double r0, double* expo, int* np, double* vectors, char* tipo) {
    int orbital;
    double sum, partial, occ, pi;
    extern double orbital_gto_imp(int , int , double , double , double* , int* , double* );

    pi = 4.f*atan(1.f);
    if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
      occ = 2.f;
    else
      occ = 1.f;

    sum = 0.f;
    for (orbital = 0; orbital < elec; orbital++) {
      partial = orbital_gto_imp(nt, orbital, r, r0, expo, np, vectors);
      sum = sum + occ*partial*partial;
    }

    return(sum/(4.f*pi));
 }

  double der_rho_radial_gto_imp(int nt, int elec, double r, double r0, double* expo, int* np, double* vectors, char* tipo) {
    int orbital;
    double sum, partial, occ, pi;
    extern double orbital_gto_imp(int , int , double , double , double* , int* , double* );
    extern double der_orbital_gto_imp(int , int , double , double , double* , int* , double* );

    pi = 4.f*atan(1.f);
    if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
      occ = 2.f;
    else
      occ = 1.f;

    sum = 0.f;
    for (orbital = 0; orbital < elec; orbital++) {
      partial = 2.f*orbital_gto_imp(nt, orbital, r, r0, expo, np, vectors);
      partial = partial*der_orbital_gto_imp(nt, orbital, r, r0, expo, np, vectors);
      sum = sum + occ*partial;
    }

    return(sum/(4.f*pi));
 }  

  double sec_der_rho_radial_gto_imp(int nt, int elec, double r, double r0, double* expo, int* np, double* vectors, char* tipo) {
    int orbital;
    double sum, partial, partial1, partial2, occ, pi;
    extern double orbital_gto_imp(int , int , double , double , double* , int* , double* );
    extern double der_orbital_gto_imp(int , int , double , double , double* , int* , double* );
    extern double sec_der_orbital_gto_imp(int , int , double , double , double* , int* , double* );

    pi = 4.f*atan(1.f);
    if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
      occ = 2.f;
    else
      occ = 1.f;

    sum = 0.f;
    for (orbital = 0; orbital < elec; orbital++) {
      partial1 = orbital_gto_imp(nt, orbital, r, r0, expo, np, vectors);
      partial1 = partial1*sec_der_orbital_gto_imp(nt, orbital, r, r0, expo, np, vectors);
 
      partial2 = der_orbital_gto_imp(nt, orbital, r, r0, expo, np, vectors);
      partial2 = partial2*der_orbital_gto_imp(nt, orbital, r, r0, expo, np, vectors);
      
      partial = 2.f*(partial1 + partial2);
      
      sum = sum + occ*partial;
    }

    return(sum/(4.f*pi));
 }

