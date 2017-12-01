#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>



 double msto_finite_int(char   *using_gamma,
                        int     mu, 
                        double  r, 
                        double  Rc, 
                        int    *np, 
                        double  gamma_couple, 
                        double *zeta, 
                        double *N_minus)

 {
double total;

 if (strcmp(using_gamma,"YES") == 0) {//label 1
 
   //total= N_minus[mu]*(Rc - r*gamma_couple)*exp(-zeta[mu]*r);
   if(r == 0){
     if (np[mu] == 1)
       total = N_minus[mu]*Rc;
     else
       total = 0.f;
   }
   else
   {
     if (np[mu] == 1)
       total = 1.f;
     else
       total = pow(r, (double) np[mu] - 1);
   
     total = total*N_minus[mu]*(Rc - r*gamma_couple);
     total = total*exp(-r*zeta[mu]);
   }

  }//label 1
   else {//label 2

        if(r == 0){
          if (np[mu] == 1)
            total = N_minus[mu];
          else
            total = 0.f;
        }
        else
        { 
          if (np[mu] == 1)
            total = 1.f;
          else
            total = pow(r, (double) np[mu] - 1);
          
          total = total*N_minus[mu];
          total = total*exp(-r*zeta[mu]);
        }

       }//label 2

  return (total);
 }


 double msto_finite_ext(int     mu, 
                        double  r, 
                        double  Rc, 
                        int    *mang, 
                        double  gamma_couple, 
                        double *alfa, 
                        double *N_plus)

 {
 double total;
 

    total= N_plus[mu]*exp(-alfa[mu]*r);
 
    total = total*pow(r,(double) -(mang[mu] + 1));

  return (total);
 }

//Derivada_MSTO_parte_interna "dfmu/dr"
 double der_msto_finite_int(char   *using_gamma,
                            int     mu, 
                            double  r, 
                            double  Rc, 
                            int    *np, 
                            double  gamma_couple, 
                            double *zeta, 
                            double *N_minus)
 {
 extern double msto_finite_int(char   *using_gamma,
                               int     mu, 
                               double  r, 
                               double  Rc, 
                               int    *np, 
                               double  gamma_couple, 
                               double *zeta, 
                               double *N_minus);
 double total;

 if (strcmp(using_gamma,"YES") == 0) {//label 3

	if(r == 0){

	if(np[mu] == 1)  
                        total = -N_minus[mu]*(gamma_couple + zeta[mu]*Rc);
          else
	    if(np[mu] == 2) 
                            total = N_minus[mu]*Rc;
                            
             else
                 total = 0.f;
	          }
	else
 	{
        if(np[mu] == 1){
                    total = N_minus[mu]*(gamma_couple + Rc*zeta[mu] - r*gamma_couple*zeta[mu]);
                        total = -total*exp(-zeta[mu]*r);
                       }
                      
           else
             if(np[mu] == 2)
                      {
                       total = N_minus[mu]*(Rc - r*Rc*zeta[mu] + r*gamma_couple*(r*zeta[mu] -2.f));
                         total = total*exp(-zeta[mu]*r);
                      }
                   else
                        {
                         total =(double) N_minus[mu]*(Rc*(np[mu] - 1.f - r*zeta[mu]) + r*gamma_couple*(r*zeta[mu] - np[mu]));
                            total = total*exp(-zeta[mu]*r)*pow(r, (double) np[mu] - 2);

                        }
         }

  }//label 3
   else {//label 4

        if(r == 0){

        if(np[mu] == 1)
                        total = -N_minus[mu]*zeta[mu];
          else
            if(np[mu] == 2)
                            total = N_minus[mu];
             else
                 total = 0.f;
                  }
        else
        {
        if(np[mu] == 1){
                    total = N_minus[mu]*zeta[mu];
                        total = -total*exp(-zeta[mu]*r);
                       }

           else
             if(np[mu] == 2)
                      {
                       total = N_minus[mu]*(r*zeta[mu] - 1.f);
                         total = -total*exp(-zeta[mu]*r);
                      }
                   else
                        {
                         total =(double) N_minus[mu]*(1 + r*zeta[mu] - np[mu]);
                            total = -total*exp(-zeta[mu]*r)*pow(r, (double) np[mu] - 2);

                        }
         }



       }//label 4


 return (total);

}
//Segunda derivada del msto parte interna  "d²fmu/dr²"
 double sec_der_msto_finite_int(char   *using_gamma,
                                int     mu,
                                double  r,
                                double  Rc,
                                int    *np,
                                double  gamma_couple,
                                double *zeta,
                                double *N_minus)
 {
 extern double msto_finite_int(char   *using_gamma,
                               int     mu,
                               double  r,
                               double  Rc,
                               int    *np,
                               double  gamma_couple,
                               double *zeta,
                               double *N_minus);
 double total;

 if (strcmp(using_gamma,"YES") == 0) {//label 3

   printf("\nNO HAY EXPRESIONES\n");
  }//label 3
   else {//label 4

        if(r == 0.f){

        if (np[mu] == 1)
           total = N_minus[mu]*pow(zeta[mu],2.f);
          else
            if(np[mu] == 2)
              total = -2.f*N_minus[mu]*zeta[mu];
             else
               if(np[mu] == 3)
                 total = 2.f*N_minus[mu];
                else 
                 total = 0.f;
        }
        else
        {
         total = N_minus[mu]*exp(-zeta[mu]*r);
         if(np[mu] == 1)
             total = total*pow(zeta[mu],2.f);
            else
              if(np[mu] == 2)
                 total = total*(-2.f*zeta[mu] + r*pow(zeta[mu],2.f));
                else
                  if(np[mu] == 3)
                     total = total*(pow(r,2.f)*pow(zeta[mu],2.f) + 2.f - 4.f*r*zeta[mu]);
                    else
                         {
                          total =total*(pow(r,(double) (np[mu] - 1))*pow(zeta[mu],2.f) 
                                  + (double) (np[mu] - 1)*pow(r, (double) (np[mu] - 3))*(np[mu] - 2.f - 2.f*r*zeta[mu]));

                         }
         }



       }//label 4


 return (total);

}

//Derivada_MSTO_parte_externa "dfmu/dr"
 double der_msto_finite_ext(int     mu, 
                            double  r, 
                            double  Rc, 
                            int    *mang, 
                            double  gamma_couple, 
                            double *alfa, 
                            double *N_plus)

 {
double total;


total = (double) alfa[mu]*r + mang[mu] + 1.f;
total = total*pow(r, (double) -(mang[mu] + 2));
total = total*exp(-alfa[mu]*r);
total = -total*N_plus[mu];

return (total);


}
//Segunda derivada del msto parte externa "d²fmu/dr²"
 double sec_der_msto_finite_ext(int     mu,
                                double  r,
                                double  Rc,
                                int    *mang,        //lmu
                                double  gamma_couple,
                                double *alfa,
                                double *N_plus)

{
double total;

// total = N_plus[mu]*pow(r, (double) -(mang[mu] + 1.f))*exp(-alfa[mu]*r);
// total =   (double) ((mang[mu] + 1)*(mang[mu] + 2) + (2*(mang[mu] + 1) + alfa[mu]*r)*alfa[mu]*r)*pow(r,-2.f)*total;
   total = N_plus[mu]*pow(r, (double) -(mang[mu] + 3.f))*exp(-alfa[mu]*r);
   total = total*((double) 2.f + pow((double) mang[mu],2.f) + 2.f*r*alfa[mu] + pow((double) r*alfa[mu],2.f) + mang[mu]*(3.f + 2.f*r*alfa[mu]));


return (total);

}

//Orbital parte interna
 double finite_orbital_int(char   *using_gamma,
                           int     nt, 
                           int     orbital, 
                           double  r, 
                           double  Rc, 
                           double *zeta,
                           int    *np, 
                           double *vectors, 
                           double  gamma_couple, 
                           double *N_minus)
 {
  int i;
  double suma;
  extern double msto_finite_int(char   *using_gamma,
                                int     mu, 
                                double  r, 
                                double  Rc, 
                                int    *np, 
                                double  gamma_couple, 
                                double *zeta, 
                                double *N_minus);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*msto_finite_int(using_gamma,
                                                          i, 
                                                          r, 
                                                          Rc, 
                                                          np, 
                                                          gamma_couple, 
                                                          zeta, 
                                                          N_minus);
  }

  return (suma);
  }
//Orbital parte externa
 double finite_orbital_ext(int     nt, 
                           int     orbital, 
                           double  r, 
                           double  Rc, 
                           double *alfa,
                           int    *mang, 
                           double *vectors, 
                           double  gamma_couple, 
                           double *N_plus)
 {
  int i;
  double suma;
  extern double msto_finite_ext(int     mu, 
                                double  r, 
                                double  Rc, 
                                int    *mang, 
                                double  gamma_couple, 
                                double *alfa, 
                                double *N_plus);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*msto_finite_ext(i, 
                                                          r, 
                                                          Rc, 
                                                          mang, 
                                                          gamma_couple, 
                                                          alfa, 
                                                          N_plus);
  }

  return (suma);
  }
//Derivada del orbital parte interna
 double der_finite_orbital_int(char   *using_gamma,
                               int     nt, 
                               int     orbital, 
                               double  r, 
                               double  Rc, 
                               double *zeta,
                               int    *np, 
                               double *vectors, 
                               double  gamma_couple, 
                               double *N_minus)
 {
  int i;
  double suma;
  extern double der_msto_finite_int(char   *using_gamma,
                                    int     mu, 
                                    double  r, 
                                    double  Rc, 
                                    int    *np, 
                                    double  gamma_couple, 
                                    double *zeta, 
                                    double *N_minus);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*der_msto_finite_int(using_gamma,
                                                              i, 
                                                              r, 
                                                              Rc, 
                                                              np, 
                                                              gamma_couple, 
                                                              zeta, 
                                                              N_minus);
  }

  return (suma);
  }
//Derivada del orbital parte externa
 double der_finite_orbital_ext(int     nt, 
                               int     orbital, 
                               double  r, 
                               double  Rc, 
                               double *alfa,
                               int    *mang, 
                               double *vectors, 
                               double  gamma_couple, 
                               double *N_plus)
 {
  int i;
  double suma;
  extern double der_msto_finite_ext(int     mu, 
                                    double  r, 
                                    double  Rc, 
                                    int    *mang, 
                                    double  gamma_couple, 
                                    double *alfa, 
                                    double *N_plus);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*der_msto_finite_ext(i,       
                                                              r,       
                                                              Rc,       
                                                              mang,       
                                                              gamma_couple,       
                                                              alfa,       
                                                              N_plus);
  }

  return (suma);
  }
//////////////////////////////////////////////////////////////////////////////
//Segunda derivada del orbital parte interna
 double sec_der_finite_orbital_int(char   *using_gamma,
                                   int     nt,
                                   int     orbital,
                                   double  r,
                                   double  Rc,
                                   double *zeta,
                                   int    *np,
                                   double *vectors,
                                   double  gamma_couple,
                                   double *N_minus)
 {
  int i;
  double suma;
  extern double sec_der_msto_finite_int(char   *using_gamma,
                                        int     mu,
                                        double  r,
                                        double  Rc,
                                        int    *np,
                                        double  gamma_couple,
                                        double *zeta,
                                        double *N_minus);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*sec_der_msto_finite_int(using_gamma,
                                                                  i,
                                                                  r,
                                                                  Rc,
                                                                  np,
                                                                  gamma_couple,
                                                                  zeta,
                                                                  N_minus);
  }

  return (suma);
  }
//////////////////////////////////////////////////////////////////////////////
//Segunda derivada del orbital parte externa
 double sec_der_finite_orbital_ext(int     nt,
                                   int     orbital,
                                   double  r,
                                   double  Rc,
                                   double *alfa,
                                   int    *mang,
                                   double *vectors,
                                   double  gamma_couple,
                                   double *N_plus)

 {
  int i;
  double suma;
  extern  double sec_der_msto_finite_ext(int     mu,
                                         double  r,
                                         double  Rc,
                                         int    *mang,        //lmu
                                         double  gamma_couple,
                                         double *alfa,
                                         double *N_plus);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*sec_der_msto_finite_ext(i,
                                                                  r,
                                                                  Rc,
                                                                  mang,        //lmu
                                                                  gamma_couple,
                                                                  alfa,
                                                                  N_plus);

  }

  return (suma);
  }

//////////////////////////////////////////////////////////////////////////////

//Densidad orbital radial parte interna
 double rho_orbital_radial_finite_int(char   *using_gamma,
                              int     nt, 
                              int     elec, 
                              double  r, 
                              double  Rc, 
                              double *zeta,
                              int    *np, 
                              double *vectors, 
                              char   *tipo, 
                              double  gamma_couple, 
                              double *N_minus,
                              int selected_orb)
 {
  int orbital;
  double sum, partial, occ, pi;
  extern double finite_orbital_int(char   *using_gamma,
                                   int     nt, 
                                   int     orbital, 
                                   double  r, 
                                   double  Rc, 
                                   double *zeta,
                                   int    *np, 
                                   double *vectors, 
                                   double  gamma_couple, 
                                   double *N_minus);


  pi = 4.f*atan(1.f);
  sum = 0.f;
  partial = finite_orbital_int(using_gamma,
                                       nt, 
                                       selected_orb, 
                                       r, 
                                       Rc, 
                                       zeta, 
                                       np, 
                                       vectors, 
                                       gamma_couple, 
                                       N_minus);
  sum = partial*partial;
  return(sum/(4.f*pi));
 }
//Densidad orbital radial parte externa
 double rho_orbital_radial_finite_ext(int     nt, 
                              int     elec, 
                              double  r, 
                              double  Rc, 
                              double *alfa,
                              int    *mang, 
                              double *vectors, 
                              char   *tipo, 
                              double  gamma_couple, 
                              double *N_plus,
                              int selected_orb)
 {
  int orbital;
  double sum, partial, occ, pi;
  extern double finite_orbital_ext(int     nt, 
                                   int     orbital, 
                                   double  r, 
                                   double  Rc, 
                                   double *alfa,
                                   int    *mang, 
                                   double *vectors, 
                                   double  gamma_couple, 
                                   double *N_plus);




  pi = 4.f*atan(1.f);
  sum = 0.f;
  partial = finite_orbital_ext(nt, 
                               selected_orb, 
                               r, 
                                       Rc, 
                                       alfa, 
                                       mang, 
                                       vectors, 
                                       gamma_couple, 
                                       N_plus);

  sum = partial*partial;

  return(sum/(4.f*pi));
 }
//////////////////////////////////////////////////////////////////////////////

//Densidad radial parte interna
 double rho_radial_finite_int(char   *using_gamma,
                              int     nt, 
                              int     elec, 
                              double  r, 
                              double  Rc, 
                              double *zeta,
                              int    *np, 
                              double *vectors, 
                              char   *tipo, 
                              double  gamma_couple, 
                              double *N_minus)
 {
  int orbital;
  double sum, partial, occ, pi;
  extern double finite_orbital_int(char   *using_gamma,
                                   int     nt, 
                                   int     orbital, 
                                   double  r, 
                                   double  Rc, 
                                   double *zeta,
                                   int    *np, 
                                   double *vectors, 
                                   double  gamma_couple, 
                                   double *N_minus);


  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {

          partial = finite_orbital_int(using_gamma,
                                       nt, 
                                       orbital, 
                                       r, 
                                       Rc, 
                                       zeta, 
                                       np, 
                                       vectors, 
                                       gamma_couple, 
                                       N_minus);

    sum = sum + occ*partial*partial;
  }

  return(sum/(4.f*pi));
 }
//Densidad radial parte externa
 double rho_radial_finite_ext(int     nt, 
                              int     elec, 
                              double  r, 
                              double  Rc, 
                              double *alfa,
                              int    *mang, 
                              double *vectors, 
                              char   *tipo, 
                              double  gamma_couple, 
                              double *N_plus)
 {
  int orbital;
  double sum, partial, occ, pi;
  extern double finite_orbital_ext(int     nt, 
                                   int     orbital, 
                                   double  r, 
                                   double  Rc, 
                                   double *alfa,
                                   int    *mang, 
                                   double *vectors, 
                                   double  gamma_couple, 
                                   double *N_plus);




  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {

          partial = finite_orbital_ext(nt, 
                                       orbital, 
                                       r, 
                                       Rc, 
                                       alfa, 
                                       mang, 
                                       vectors, 
                                       gamma_couple, 
                                       N_plus);

    sum = sum + occ*partial*partial;
  }

  return(sum/(4.f*pi));
 }
//Derivada de la densidad radial parte interna
 double der_rho_radial_finite_int(char   *using_gamma,
                                  int     nt, 
                                  int     elec, 
                                  double  r, 
                                  double  Rc, 
                                  double *zeta,
                                  int    *np, 
                                  double *vectors, 
                                  char   *tipo, 
                                  double  gamma_couple, 
                                  double *N_minus)
 {
  int orbital;
  double sum, partial, occ, pi;
  extern double der_finite_orbital_int(char   *using_gamma,
                                       int     nt, 
                                       int     orbital, 
                                       double  r, 
                                       double  Rc, 
                                       double *zeta,
                                       int    *np, 
                                       double *vectors, 
                                       double  gamma_couple, 
                                       double *N_minus);

  extern double finite_orbital_int(char   *using_gamma,
                                   int     nt, 
                                   int     orbital, 
                                   double  r, 
                                   double  Rc, 
                                   double *zeta,
                                   int    *np, 
                                   double *vectors, 
                                   double  gamma_couple, 
                                   double *N_minus);


  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {
         partial = 2.f*finite_orbital_int(using_gamma,
                                          nt, 
                                          orbital, 
                                          r, 
                                          Rc, 
                                          zeta, 
                                          np, 
                                          vectors, 
                                          gamma_couple, 
                                          N_minus);

          partial = partial*der_finite_orbital_int(using_gamma,
                                                   nt, 
                                                   orbital, 
                                                   r, 
                                                   Rc, 
                                                   zeta, 
                                                   np, 
                                                   vectors, 
                                                   gamma_couple, 
                                                   N_minus);

    sum = sum + occ*partial;
  }

  return(sum/(4.f*pi));
 }
//////////////////////////////////////////////////////////////////////////////
//Segunda derivada de la densidad radial parte interna
 double sec_der_rho_radial_finite_int(char   *using_gamma,
                                      int     nt, 
                                      int     elec, 
                                      double  r, 
                                      double  Rc, 
                                      double *zeta,
                                      int    *np, 
                                      double *vectors, 
                                      char   *tipo, 
                                      double  gamma_couple, 
                                      double *N_minus)
 {
  int orbital;
  double sum, partial1, partial2, partial, occ, pi;
  extern  double sec_der_finite_orbital_int(char   *using_gamma,
                                            int     nt,
                                            int     orbital,
                                            double  r,
                                            double  Rc,
                                            double *zeta,
                                            int    *np,
                                            double *vectors,
                                            double  gamma_couple,
                                            double *N_minus); 
 
  extern double der_finite_orbital_int(char   *using_gamma,
                                       int     nt, 
                                       int     orbital, 
                                       double  r, 
                                       double  Rc, 
                                       double *zeta,
                                       int    *np, 
                                       double *vectors, 
                                       double  gamma_couple, 
                                       double *N_minus);

  extern double finite_orbital_int(char   *using_gamma,
                                   int     nt, 
                                   int     orbital, 
                                   double  r, 
                                   double  Rc, 
                                   double *zeta,
                                   int    *np, 
                                   double *vectors, 
                                   double  gamma_couple, 
                                   double *N_minus);


  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {
         partial1 = der_finite_orbital_int(using_gamma,
                                           nt, 
                                           orbital, 
                                           r, 
                                           Rc, 
                                           zeta, 
                                           np, 
                                           vectors, 
                                           gamma_couple, 
                                           N_minus);
         partial1 = pow(partial1,2.f);
         
         partial2 = sec_der_finite_orbital_int(using_gamma,
                                               nt,
                                               orbital,
                                               r,
                                               Rc,
                                               zeta,
                                               np,
                                               vectors,
                                               gamma_couple,
                                               N_minus);
  
         partial2 = partial2*finite_orbital_int(using_gamma,
                                                nt,
                                                orbital,
                                                r,
                                                Rc,
                                                zeta,
                                                np,
                                                vectors,
                                                gamma_couple,
                                                N_minus);

         partial = 2.f*(partial1 + partial2);

    sum = sum + occ*partial;
  }

  return(sum/(4.f*pi));
 }
//////////////////////////////////////////////////////////////////////////////
//Derivada de la densidad radial parte externa
 double der_rho_radial_finite_ext(int     nt, 
                                  int     elec, 
                                  double  r, 
                                  double  Rc, 
                                  double *alfa,
                                  int    *mang, 
                                  double *vectors, 
                                  char   *tipo, 
                                  double  gamma_couple, 
                                  double *N_plus)
 {
  int orbital;
  double sum, partial, occ, pi;
  extern double der_finite_orbital_ext(int     nt, 
                                       int     orbital, 
                                       double  r, 
                                       double  Rc, 
                                       double *alfa,
                                       int    *mang, 
                                       double *vectors, 
                                       double  gamma_couple, 
                                       double *N_plus);

  extern double finite_orbital_ext(int     nt, 
                                   int     orbital, 
                                   double  r, 
                                   double  Rc, 
                                   double *alfa,
                                   int    *mang, 
                                   double *vectors, 
                                   double  gamma_couple, 
                                   double *N_plus);




  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {
         partial = 2.f*finite_orbital_ext(nt, 
                                          orbital, 
                                          r, 
                                          Rc, 
                                          alfa, 
                                          mang, 
                                          vectors, 
                                          gamma_couple, 
                                          N_plus);

          partial = partial*der_finite_orbital_ext(nt, 
                                                   orbital, 
                                                   r, 
                                                   Rc, 
                                                   alfa, 
                                                   mang, 
                                                   vectors, 
                                                   gamma_couple, 
                                                   N_plus);

    sum = sum + occ*partial;
  }

  return(sum/(4.f*pi));
 }
//////////////////////////////////////////////////////////////////////////////
//Segunda derivada de la densidad radial parte externa
 double sec_der_rho_radial_finite_ext(int     nt, 
                                      int     elec, 
                                      double  r, 
                                      double  Rc, 
                                      double *alfa,
                                      int    *mang, 
                                      double *vectors, 
                                      char   *tipo, 
                                      double  gamma_couple, 
                                      double *N_plus)
 {
  int orbital;
  double sum, partial, partial1, partial2, occ, pi;
  extern double sec_der_finite_orbital_ext(int     nt,
                                           int     orbital,
                                           double  r,
                                           double  Rc,
                                           double *alfa,
                                           int    *mang,
                                           double *vectors,
                                           double  gamma_couple,
                                           double *N_plus);

  extern double der_finite_orbital_ext(int     nt, 
                                       int     orbital, 
                                       double  r, 
                                       double  Rc, 
                                       double *alfa,
                                       int    *mang, 
                                       double *vectors, 
                                       double  gamma_couple, 
                                       double *N_plus);

  extern double finite_orbital_ext(int     nt, 
                                   int     orbital, 
                                   double  r, 
                                   double  Rc, 
                                   double *alfa,
                                   int    *mang, 
                                   double *vectors, 
                                   double  gamma_couple, 
                                   double *N_plus);

  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {
         partial1 = der_finite_orbital_ext(nt,
                                           orbital,
                                           r,
                                           Rc,
                                           alfa,
                                           mang,
                                           vectors,
                                           gamma_couple,
                                           N_plus);
  
         partial1 = pow(partial1,2.f);

         partial2 = sec_der_finite_orbital_ext(nt,
                                               orbital,
                                               r,
                                               Rc,
                                               alfa,
                                               mang,
                                               vectors,
                                               gamma_couple,
                                               N_plus);

         partial2 = partial2*finite_orbital_ext(nt,
                                                orbital,
                                                r,
                                                Rc,
                                                alfa,
                                                mang,
                                                vectors,
                                                gamma_couple,
                                                N_plus);
 
         
         partial = 2.f*(partial1 + partial2); 


    sum = sum + occ*partial;
  }

  return(sum/(4.f*pi));
 }

