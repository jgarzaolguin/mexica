#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>



 void  print_spherical_density_confined(int     compara,
                                        int     nt,
                                        double  z,
                                        int     elecalfa,
                                        int     elecbeta,
                                        double  Rc,
                                        double *expo,
                                        int    *np,
                                        int    *mang,
                                        double *vectsfinalfa,
                                        double *vectsfinbeta,
                                        char   *tipo,
                                        double *arreglo_factorial,
                                        double *arreglo_inv_factorial)

 {//label 1

 extern double rho_radial_confined(int     nt,
                                  int      elec,
                                  double   r,
                                  double   Rc,
                                  double  *expo,
                                  int     *np,
                                  double  *vectors,
                                  char    *tipo,
                                  double  *arreglo_factorial,
                                  double  *arreglo_inv_factorial);

 extern double der_rho_radial_confined(int     nt,
                                       int     elec,
                                       double  r,
                                       double  Rc,
                                       double *expo,
                                       int    *np,
                                       double *vectors,
                                       char   *tipo,
                                       double *arreglo_factorial,
                                       double *arreglo_inv_factorial);


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
      sprintf(name, "confined_rhf_rho_drho_+drho_rdf_divdrhorho_Rc_%.4f", Rc);
   write_out = fopen(name, "w");

   do{

      rho =  rho_radial_confined(nt,
                                 elecalfa,
                                 radius,
                                 Rc,
                                 expo,
                                 np,
                                 vectsfinalfa,
                                 tipo,
                                 arreglo_factorial,
                                 arreglo_inv_factorial);
      
      der_rho =  der_rho_radial_confined(nt,
                                         elecalfa,
                                         radius,
                                         Rc,
                                         expo,
                                         np,
                                         vectsfinalfa,
                                         tipo,
                                         arreglo_factorial,
                                         arreglo_inv_factorial);
      

                  rdf     = 4.f*numero_pi*radius*radius*rho;
                  div_rho = der_rho/rho;

              /*    fprintf(write_out,"%5.4f  %24.16f  %24.16f  %24.16f  %24.16f %24.16f %24.16f\n", radius,
                                                                                                   rho,
                                                                                                   der_rho,
                                                                                                   -1.f*der_rho,
                                                                                                   rdf,
                                                                                                   div_rho,
                                                                                                   -1.f*div_rho/(2.f*z));*/
//Agregado Adrian

                  fprintf(write_out,"%5.4f  %24.16f  %24.16f\n", radius,rho,der_rho);
                                                                                              







//mrb steps            }//End for n_points
          radius = radius + 0.0001;
          if (rho < rho_ini) rho_ini = rho;
          else  rho = 0.f;

        } while(radius < Rc || fabs(rho) > bound_rho);


   }//label 2
    else {//label 3

            sprintf(name, "confined_uhf_rho_drho_+drho_rdf_divdrhorho_Rc_%f", Rc);

            write_out = fopen(name, "w");


        do{


           rho =  rho_radial_confined(nt,
                                      elecalfa,
                                      radius,
                                      Rc,
                                      expo,
                                      np,
                                      vectsfinalfa,
                                      tipo,
                                      arreglo_factorial,
                                      arreglo_inv_factorial);

           rho =  rho + rho_radial_confined(nt,
                                            elecbeta,
                                            radius,
                                            Rc,
                                            expo,
                                            np,
                                            vectsfinbeta,
                                            tipo,
                                            arreglo_factorial,
                                            arreglo_inv_factorial);

      
           der_rho =  der_rho_radial_confined(nt,
                                              elecalfa,
                                              radius,
                                              Rc,
                                              expo,
                                              np,
                                              vectsfinalfa,
                                              tipo,
                                              arreglo_factorial,
                                              arreglo_inv_factorial);

           der_rho =  der_rho + der_rho_radial_confined(nt,
                                                        elecbeta,
                                                        radius,
                                                        Rc,
                                                        expo,
                                                        np,
                                                        vectsfinbeta,
                                                        tipo,
                                                        arreglo_factorial,
                                                        arreglo_inv_factorial);

      


                  rdf     = 4.f*numero_pi*radius*radius*rho;
                  div_rho = der_rho/rho;

                  fprintf(write_out,"%5.4f  %32.16f  %32.16f  %32.16f  %32.16f %32.16f %32.16f\n", radius,rho,
                                                                                                   der_rho,
                                                                                                   -1.f*der_rho,
                                                                                                   rdf,
                                                                                                   div_rho,
                                                                                                   -1.f*div_rho/(2.f*z));



//mrb steps            }//End for n_points

          radius = radius + 0.0001;
          if (rho < rho_ini) rho_ini = rho;
          else rho = 0.f;

        } while(radius < Rc || fabs(rho) > bound_rho);


        } //label 3


       fclose(write_out);



 }//label 1
