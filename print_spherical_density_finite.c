#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>



 void  print_spherical_density_finite(int     compara,
                                      char   *using_gamma,
                                      int     nt,
                                      double  z,
                                      int     elecalfa,
                                      int     elecbeta,
                                      double  Rc,
                                      double *expo,
                                      double *zetas,
                                      int    *np,
                                      int    *mang,
                                      double *vectsfinalfa,
                                      double *vectsfinbeta,
                                      char   *tipo,
                                      double  gamma_couple,
                                      double *NC_minus,
                                      double *NC_plus)


 {//label 1


 extern double  rho_radial_finite_int(char   *using_gamma,
                                      int     nt,
                                      int     elec,
                                      double  r,
                                      double  Rc,
                                      double *zeta,
                                      int    *np,
                                      double *vectors,
                                      char   *tipo,
                                      double  gamma_couple,
                                      double *N_minus);

 extern double  rho_radial_finite_ext(int     nt,
                                      int     elec,
                                      double  r,
                                      double  Rc,
                                      double *alfa,
                                      int    *mang,
                                      double *vectors,
                                      char   *tipo,
                                      double  gamma_couple,
                                      double *N_plus);


 extern double  der_rho_radial_finite_int(char   *using_gamma,
                                          int     nt,
                                          int     elec,
                                          double  r,
                                          double  Rc,
                                          double *zeta,
                                          int    *np,
                                          double *vectors,
                                          char   *tipo,
                                          double  gamma_couple,
                                          double *N_minus);

 extern double  der_rho_radial_finite_ext(int     nt,
                                          int     elec,
                                          double  r,
                                          double  Rc,
                                          double *alfa,
                                          int    *mang,
                                          double *vectors,
                                          char   *tipo,
                                          double  gamma_couple,
                                          double *N_plus);



 char name[200];
 FILE *write_out;
 double radius,
        bound_rho,
        rho_ini,
        rho,
        der_rho,
        rdf,
        div_rho,
        numero_pi;

 double doble1;

 doble1 = (double)1;
 numero_pi = atan(doble1)*(double)4;

        radius    = 0.f;
        bound_rho = 1e-10;
        rho_ini   = 1e12;

 if (compara == 0) {//label 2
   if (strcmp(using_gamma,"YES")==0)
      sprintf(name, "rhf_rho_drho_+drho_rdf_divdrhorho_Rc_%.4f_%.4f", Rc, gamma_couple);
   else
      sprintf(name, "rhf_rho_drho_+drho_rdf_divdrhorho_Rc_%.4f_no_gamma", Rc);
   write_out = fopen(name, "w");
   do{
     if(radius < Rc) {
       rho     = rho_radial_finite_int(using_gamma,
                                                 nt,
                                                 elecalfa,
                                                 radius,
                                                 Rc,
                                                 expo,
                                                 np,
                                                 vectsfinalfa,
                                                 tipo,
                                                 gamma_couple,
                                                 NC_minus);

                  der_rho = der_rho_radial_finite_int(using_gamma,
                                                      nt,
                                                      elecalfa,
                                                      radius,
                                                      Rc,
                                                      expo,
                                                      np,
                                                      vectsfinalfa,
                                                      tipo,
                                                      gamma_couple,
                                                      NC_minus);


                  rdf     = 4.f*numero_pi*radius*radius*rho;
                  div_rho = der_rho/rho;

                  fprintf(write_out,"%5.4f  %24.16f  %24.16f  %24.16f  %24.16f %24.16f %24.16f\n", radius,
                                                                                                   rho,
                                                                                                   der_rho,
                                                                                                   -1.f*der_rho,
                                                                                                   rdf,
                                                                                                   div_rho,
                                                                                                   -1.f*div_rho/(2.f*z));


               }
               if(radius >= Rc)
               {
                 rho     = rho_radial_finite_ext(nt,
                                                 elecalfa,
                                                 radius,
                                                 Rc,
                                                 zetas,
                                                 mang,
                                                 vectsfinalfa,
                                                 tipo,
                                                 gamma_couple,
                                                 NC_plus);

                  der_rho = der_rho_radial_finite_ext(nt,
                                                      elecalfa,
                                                      radius,
                                                      Rc,
                                                      zetas,
                                                      mang,
                                                      vectsfinalfa,
                                                      tipo,
                                                      gamma_couple,
                                                      NC_plus);

                  rdf     = 4.f*numero_pi*radius*radius*rho;
                  div_rho = der_rho/rho;

                   fprintf(write_out,"%5.4f  %24.16f  %24.16f %24.16f  %24.16f %24.16f %24.16f\n", radius,
                                                                                                   rho,
                                                                                                   der_rho,
                                                                                                   -1.f*der_rho,
                                                                                                   rdf,
                                                                                                   div_rho,
                                                                                                   -1.f*div_rho/(2.f*z));


               }


//mrb steps            }//End for n_points
          radius = radius + 0.0001;
          if (rho < rho_ini) rho_ini = rho;
          else  rho = 0.f;

        }while(fabs(rho) > bound_rho);


   }//label 2
    else {//label 3

      if (strcmp(using_gamma,"YES")==0)
            sprintf(name, "uhf_rho_drho_+drho_rdf_divdrhorho_Rc_%f_%f", Rc, gamma_couple);
         else
            sprintf(name, "uhf_rho_drho_+drho_rdf_divdrhorho_Rc_%f_no_gama", Rc);

            write_out = fopen(name, "w");


        do{

              if(radius < Rc)
               {
                 rho     = rho_radial_finite_int(using_gamma,
                                                 nt,
                                                 elecalfa,
                                                 radius,
                                                 Rc,
                                                 expo,
                                                 np,
                                                 vectsfinalfa,
                                                 tipo,
                                                 gamma_couple,
                                                 NC_minus);

                 rho         = rho + rho_radial_finite_int(using_gamma,
                                                           nt,
                                                           elecbeta,
                                                           radius,
                                                           Rc,
                                                           expo,
                                                           np,
                                                           vectsfinbeta,
                                                           tipo,
                                                           gamma_couple,
                                                           NC_minus);

                  der_rho = der_rho_radial_finite_int(using_gamma,
                                                      nt,
                                                      elecalfa,
                                                      radius,
                                                      Rc,
                                                      expo,
                                                      np,
                                                      vectsfinalfa,
                                                      tipo,
                                                      gamma_couple,
                                                      NC_minus);

                  der_rho = der_rho + der_rho_radial_finite_int(using_gamma,
                                                                nt,
                                                                elecbeta,
                                                                radius,
                                                                Rc,
                                                                expo,
                                                                np,
                                                                vectsfinbeta,
                                                                tipo,
                                                                gamma_couple,
                                                                NC_minus);

                  rdf     = 4.f*numero_pi*radius*radius*rho;
                  div_rho = der_rho/rho;

                  fprintf(write_out,"%5.4f  %32.16f  %32.16f  %32.16f  %32.16f %32.16f %32.16f\n", radius,
                                                                                                   rho,
                                                                                                   der_rho,
                                                                                                   -1.f*der_rho,
                                                                                                   rdf,
                                                                                                   div_rho,
                                                                                                   -1.f*div_rho/(2.f*z));



               }
               if(radius >= Rc)
               {
                 rho     = rho_radial_finite_ext(nt,
                                                 elecalfa,
                                                 radius,
                                                 Rc,
                                                 zetas,
                                                 mang,
                                                 vectsfinalfa,
                                                 tipo,
                                                 gamma_couple,
                                                 NC_plus);

                 rho         = rho + rho_radial_finite_ext(nt,
                                                           elecbeta,
                                                           radius,
                                                           Rc,
                                                           zetas,
                                                           mang,
                                                           vectsfinbeta,
                                                           tipo,
                                                           gamma_couple,
                                                           NC_plus);

                  der_rho = der_rho_radial_finite_ext(nt,
                                                      elecalfa,
                                                      radius,
                                                      Rc,
                                                      zetas,
                                                      mang,
                                                      vectsfinalfa,
                                                      tipo,
                                                      gamma_couple,
                                                      NC_plus);

                  der_rho = der_rho + der_rho_radial_finite_ext(nt,
                                                                elecbeta,
                                                                radius,
                                                                Rc,
                                                                zetas,
                                                                mang,
                                                                vectsfinbeta,
                                                                tipo,
                                                                gamma_couple,
                                                                NC_plus);


                  rdf     = 4.f*numero_pi*radius*radius*rho;
                  div_rho = der_rho/rho;

                   fprintf(write_out,"%5.4f  %32.16f  %32.16f %32.16f  %32.16f %32.16f %32.16f\n", radius,
                                                                                                   rho,
                                                                                                   der_rho,
                                                                                                   -1.f*der_rho,
                                                                                                   rdf,
                                                                                                   div_rho,
                                                                                                   -1.f*div_rho/(2.f*z));


               }
//mrb steps            }//End for n_points

          radius = radius + 0.0001;
          if (rho < rho_ini) rho_ini = rho;
          else rho = 0.f;

        }while(fabs(rho) > bound_rho);


        } //label 3


       fclose(write_out);



 }//label 1
