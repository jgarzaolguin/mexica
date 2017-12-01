#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

 int shannon_entropy_finite(int     z, 
                            char   *using_gamma,
                            int     compara,
                            int     nt,
                            int     elecalfa,
                            int     elecbeta,
                            double  Rc,
                            double *expo,
                            int    *np,
                            double *zetas,
                            int    *mang,
                            double *vectsfinalfa,
                            double *vectsfinbeta,
                            char   *tipo,
                            double *NC_minus,
                            double *NC_plus,
                            double  gamma_couple)
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



  double numero_pi;

  double doble1;

       doble1 = (double)1;
       numero_pi = atan(doble1)*(double)4;


  double  point_r,
          p1,
          p2,
          p3,
          rho_1,
          rho_2,
          rho_3;
  double  powp1, powp2, powp3, m, n, l;
  int     n_points;
  double  suma;
  double  suma1;
  int     i;
  double  sumas;
  double  rho_1s,
          rho_2s,
          rho_3s,
          ms,
          ns,
          ls;

  n_points = 3000;

   
  if (compara == 0){//label 2

           suma1 = 0.f;
           sumas = 0.f;

           for (i = 1; i <= (n_points/2 - 1); i++) {
             if (i == 1) {
                 p1 = 0.f;
                 point_r = -5.f;
                 p2 = exp(point_r)/z;               // rj                      // 
                 point_r = -5.f + 1.f/32.f;                             //
                 p3 = exp(point_r)/z;               // rj+2                    //
             }
             else {
                 point_r = -10.f + (double) (2*i - 3)/100;                           // maya tipo Froese Fischer
                 p1 = exp(point_r)/z;               // rj                      //    
                 point_r = -10.f + (double) (2*i - 2)/100;                             //
                 p2 = exp(point_r)/z;               // rj+2                    //    
                 point_r = -10.f + (double) (2*i - 1)/100;                             //
                 p3 = exp(point_r)/z;               // rj+2                    //
             }

 
             if (p1 < Rc)
                 rho_1 = rho_radial_finite_int(using_gamma,
                                               nt,
                                               elecalfa,
                                               p1,
                                               Rc,
                                               expo,
                                               np,
                                               vectsfinalfa,
                                               tipo,
                                               gamma_couple,
                                               NC_minus);
              else
                 rho_1 =  rho_radial_finite_ext(nt,
                                                elecalfa,
                                                p1,
                                                Rc,
                                                zetas,
                                                mang,
                                                vectsfinalfa,
                                                tipo,
                                                gamma_couple,
                                                NC_plus);

             if (p2 < Rc)
                 rho_2 = rho_radial_finite_int(using_gamma,
                                               nt,
                                               elecalfa,
                                               p2,
                                               Rc,
                                               expo,
                                               np,
                                               vectsfinalfa,
                                               tipo,
                                               gamma_couple,
                                               NC_minus);
             else
                 rho_2 =  rho_radial_finite_ext(nt,
                                                elecalfa,
                                                p2,
                                                Rc,
                                                zetas,
                                                mang,
                                                vectsfinalfa,
                                                tipo,
                                                gamma_couple,
                                                NC_plus);
             
              if (p3 < Rc)
                 rho_3 = rho_radial_finite_int(using_gamma,
                                               nt,
                                               elecalfa,
                                               p3,
                                               Rc,
                                               expo,
                                               np,
                                               vectsfinalfa,
                                               tipo,
                                               gamma_couple,
                                               NC_minus);
               else
                 rho_3 =  rho_radial_finite_ext(nt,
                                                elecalfa,
                                                p3,
                                                Rc,
                                                zetas,
                                                mang,
                                                vectsfinalfa,
                                                tipo,
                                                gamma_couple,
                                                NC_plus);


                 rho_1s = -rho_1*log(rho_1)*p1*p1;
                 rho_2s = -rho_2*log(rho_2)*p2*p2;
                 rho_3s = -rho_3*log(rho_3)*p3*p3;

                 if (rho_1s != rho_1s)
                   rho_1s = 0.f;      
                 if (rho_2s != rho_2s)
                   rho_2s = 0.f;      
                 if (rho_3s != rho_3s)
                   rho_3s = 0.f;      
               

                 rho_1 = rho_1*p1*p1;
                 rho_2 = rho_2*p2*p2;
                 rho_3 = rho_3*p3*p3;


                 powp1 = pow(p1,2.f);
                 powp2 = pow(p2,2.f);
                 powp3 = pow(p3,2.f);

                 n = ((p2 + p1)*(p3 + p2)/(p1*p2*(p3 + p2) - p2*p3*(p2 + p1)))*((powp2*rho_1 - powp1*rho_2)/(powp2 - powp1) - (powp3*rho_2 - powp2*rho_3)/(powp3 - powp2));
                 l = (1.f/(p3*(p2 + p1) - p1*(p3 + p2)))*(p3*(powp2*rho_1 - powp1*rho_2)/(p2 - p1) - p1*(powp3*rho_2 - powp2*rho_3)/(p3 - p2));
                 m = (rho_3 - n*p3 - l)/powp3;

                 ns = ((p2 + p1)*(p3 + p2)/(p1*p2*(p3 + p2) - p2*p3*(p2 + p1)))*((powp2*rho_1s - powp1*rho_2s)/(powp2 - powp1) - (powp3*rho_2s - powp2*rho_3s)/(powp3 - powp2));
                 ls = (1.f/(p3*(p2 + p1) - p1*(p3 + p2)))*(p3*(powp2*rho_1s - powp1*rho_2s)/(p2 - p1) - p1*(powp3*rho_2s - powp2*rho_3s)/(p3 - p2));
                 ms = (rho_3s - ns*p3 - ls)/powp3;

                 if (p1 < p2 && p2 < p3) {
                 
                 
                 if (rho_3 > 1e-40) {
                   suma1 = suma1 + (1.f/3.f)*m*(pow(p3,3.f) - pow(p1,3.f)) + (1.f/2.f)*n*(powp3 - powp1) + l*(p3 - p1);
                   sumas = sumas + (1.f/3.f)*ms*(pow(p3,3.f) - pow(p1,3.f)) + (1.f/2.f)*ns*(powp3 - powp1) + ls*(p3 - p1);
//                   printf("%f p1 %f p2 %f p3 %1.66f %1.66f\n", p1, p2, p3, rho_3, rho_3s);
                  } else
                     i = (n_points/2 - 1);
                 }


           }

           printf("\nRHF Shannon Entropy  & Int[rho(r),{r,0,Rc}] = %10.6f  %10.6f\n", 4.f*numero_pi*sumas, 4.f*numero_pi*suma1);



//mrb      do{
//mrb         if(radius < Rc) {
//mrb            rho = rho_radial_finite_int(using_gamma,
//mrb                                        nt,
//mrb                                        elecalfa,
//mrb                                        radius,
//mrb                                        Rc,
//mrb                                        expo,
//mrb                                        np,
//mrb                                        vectsfinalfa,
//mrb                                        tipo,
//mrb                                        gamma_couple,
//mrb                                        NC_minus);
//mrb         }
//mrb         
//mrb         if(radius >= Rc){
//mrb            rho = rho_radial_finite_ext(nt,
//mrb                                        elecalfa,
//mrb                                        radius,
//mrb                                        Rc,
//mrb                                        zetas,
//mrb                                        mang,
//mrb                                        vectsfinalfa,
//mrb                                        tipo,
//mrb                                        gamma_couple,
//mrb                                        NC_plus);
//mrb         }
//mrb         if (rho < rho_ini) {
//mrb            shannon = -rho*log(rho)*4.f*numero_pi*radius*radius;
//mrb            slope            = (shannon - shannon_i)/(radius - radius_i);
//mrb            intercept        = shannon_i - radius_i*slope;
//mrb            shannon_entropy = shannon_entropy + 
//mrb                              slope*(radius*radius - radius_i*radius_i)/2.f + 
//mrb                              intercept*(radius - radius_i);
//mrb            shannon_i = shannon;
//mrb            radius_i  = radius;
//mrb            rho_ini = rho;
//mrb            radius = 0.0001f + radius;
//mrb         }
//mrb          else {
//mrb            rho = 0.f;
//mrb            shannon = 0.f;
//mrb          }
//mrb      } while(fabs(rho) > bound_rho);

//mrb       printf("\nR_inf= %6.3f, SHANNON ENTROPY RHF %f\n", radius, shannon_entropy);

  }//label 2
   else {//label 3

           suma1 = 0.f;
           sumas = 0.f;

           for (i = 1; i <= (n_points/2 - 1); i++) {
             if (i == 1) {
                 p1 = 0.f;
                 point_r = -5.f;
                 p2 = exp(point_r)/z;               // rj                      // 
                 point_r = -5.f + 1.f/32.f;                             //
                 p3 = exp(point_r)/z;               // rj+2                    //
             }
             else {
                 point_r = -10.f + (double) (2*i - 3)/100;                           // maya tipo Froese Fischer
                 p1 = exp(point_r)/z;               // rj                      //    
                 point_r = -10.f + (double) (2*i - 2)/100;                             //
                 p2 = exp(point_r)/z;               // rj+2                    //    
                 point_r = -10.f + (double) (2*i - 1)/100;                             //
                 p3 = exp(point_r)/z;               // rj+2                    //
             }

 
             if (p1 < Rc) {
                 rho_1 = rho_radial_finite_int(using_gamma,
                                               nt,
                                               elecalfa,
                                               p1,
                                               Rc,
                                               expo,
                                               np,
                                               vectsfinalfa,
                                               tipo,
                                               gamma_couple,
                                               NC_minus);

                 rho_1 = rho_1 + rho_radial_finite_int(using_gamma,
                                                       nt,
                                                       elecbeta,
                                                       p1,
                                                       Rc,
                                                       expo,
                                                       np,
                                                       vectsfinbeta,
                                                       tipo,
                                                       gamma_couple,
                                                       NC_minus);
              }
              else {
                 rho_1 =  rho_radial_finite_ext(nt,
                                                elecalfa,
                                                p1,
                                                Rc,
                                                zetas,
                                                mang,
                                                vectsfinalfa,
                                                tipo,
                                                gamma_couple,
                                                NC_plus);

                 rho_1 =  rho_1 + rho_radial_finite_ext(nt,
                                                        elecbeta,
                                                        p1,
                                                        Rc,
                                                        zetas,
                                                        mang,
                                                        vectsfinbeta,
                                                        tipo,
                                                        gamma_couple,
                                                        NC_plus);
              }

             if (p2 < Rc) {
                 rho_2 = rho_radial_finite_int(using_gamma,
                                               nt,
                                               elecalfa,
                                               p2,
                                               Rc,
                                               expo,
                                               np,
                                               vectsfinalfa,
                                               tipo,
                                               gamma_couple,
                                               NC_minus);

                 rho_2 = rho_2 + rho_radial_finite_int(using_gamma,
                                                       nt,
                                                       elecbeta,
                                                       p2,
                                                       Rc,
                                                       expo,
                                                       np,
                                                       vectsfinbeta,
                                                       tipo,
                                                       gamma_couple,
                                                       NC_minus);
             }
             else {
                 rho_2 =  rho_radial_finite_ext(nt,
                                                elecalfa,
                                                p2,
                                                Rc,
                                                zetas,
                                                mang,
                                                vectsfinalfa,
                                                tipo,
                                                gamma_couple,
                                                NC_plus);

                 rho_2 =  rho_2 + rho_radial_finite_ext(nt,
                                                        elecbeta,
                                                        p2,
                                                        Rc,
                                                        zetas,
                                                        mang,
                                                        vectsfinbeta,
                                                        tipo,
                                                        gamma_couple,
                                                        NC_plus);
             }
             
              if (p3 < Rc) {
                 rho_3 = rho_radial_finite_int(using_gamma,
                                               nt,
                                               elecalfa,
                                               p3,
                                               Rc,
                                               expo,
                                               np,
                                               vectsfinalfa,
                                               tipo,
                                               gamma_couple,
                                               NC_minus);

                 rho_3 = rho_3 + rho_radial_finite_int(using_gamma,
                                                       nt,
                                                       elecbeta,
                                                       p3,
                                                       Rc,
                                                       expo,
                                                       np,
                                                       vectsfinbeta,
                                                       tipo,
                                                       gamma_couple,
                                                       NC_minus);
               }
               else {
                 rho_3 =  rho_radial_finite_ext(nt,
                                                elecalfa,
                                                p3,
                                                Rc,
                                                zetas,
                                                mang,
                                                vectsfinalfa,
                                                tipo,
                                                gamma_couple,
                                                NC_plus);

                 rho_3 =  rho_3 + rho_radial_finite_ext(nt,
                                                        elecbeta,
                                                        p3,
                                                        Rc,
                                                        zetas,
                                                        mang,
                                                        vectsfinbeta,
                                                        tipo,
                                                        gamma_couple,
                                                        NC_plus);
               }


                 rho_1s = -rho_1*log(rho_1)*p1*p1;
                 rho_2s = -rho_2*log(rho_2)*p2*p2;
                 rho_3s = -rho_3*log(rho_3)*p3*p3;

                 if (rho_1s != rho_1s)
                   rho_1s = 0.f;      
                 if (rho_2s != rho_2s)
                   rho_2s = 0.f;      
                 if (rho_3s != rho_3s)
                   rho_3s = 0.f;      
               

                 rho_1 = rho_1*p1*p1;
                 rho_2 = rho_2*p2*p2;
                 rho_3 = rho_3*p3*p3;


                 powp1 = pow(p1,2.f);
                 powp2 = pow(p2,2.f);
                 powp3 = pow(p3,2.f);

                 n = ((p2 + p1)*(p3 + p2)/(p1*p2*(p3 + p2) - p2*p3*(p2 + p1)))*((powp2*rho_1 - powp1*rho_2)/(powp2 - powp1) - (powp3*rho_2 - powp2*rho_3)/(powp3 - powp2));
                 l = (1.f/(p3*(p2 + p1) - p1*(p3 + p2)))*(p3*(powp2*rho_1 - powp1*rho_2)/(p2 - p1) - p1*(powp3*rho_2 - powp2*rho_3)/(p3 - p2));
                 m = (rho_3 - n*p3 - l)/powp3;

                 ns = ((p2 + p1)*(p3 + p2)/(p1*p2*(p3 + p2) - p2*p3*(p2 + p1)))*((powp2*rho_1s - powp1*rho_2s)/(powp2 - powp1) - (powp3*rho_2s - powp2*rho_3s)/(powp3 - powp2));
                 ls = (1.f/(p3*(p2 + p1) - p1*(p3 + p2)))*(p3*(powp2*rho_1s - powp1*rho_2s)/(p2 - p1) - p1*(powp3*rho_2s - powp2*rho_3s)/(p3 - p2));
                 ms = (rho_3s - ns*p3 - ls)/powp3;

                 if (p1 < p2 && p2 < p3) {
                 
                 
                 if (rho_3 > 1e-40) {
                   suma1 = suma1 + (1.f/3.f)*m*(pow(p3,3.f) - pow(p1,3.f)) + (1.f/2.f)*n*(powp3 - powp1) + l*(p3 - p1);
                   sumas = sumas + (1.f/3.f)*ms*(pow(p3,3.f) - pow(p1,3.f)) + (1.f/2.f)*ns*(powp3 - powp1) + ls*(p3 - p1);
       //            printf("%f p1 %f p2 %f p3 %1.66f %1.66f\n", p1, p2, p3, rho_3, rho_3s);
                  } else
                     i = (n_points/2 - 1);
                 }


           }

           printf("\nUHF Shannon Entropy  & Int[rho(r),{r,0,Rc}] = %10.6f  %10.6f\n", 4.f*numero_pi*sumas, 4.f*numero_pi*suma1);


//mrb      do{
//mrb        if(radius < Rc) {
//mrb          rho = rho_radial_finite_int(using_gamma,
//mrb                                      nt,
//mrb                                      elecalfa,
//mrb                                      radius,
//mrb                                      Rc,
//mrb                                      expo,
//mrb                                      np,
//mrb                                      vectsfinalfa,
//mrb                                      tipo,
//mrb                                      gamma_couple,
//mrb                                      NC_minus);
//mrb                
//mrb        rho = rho + rho_radial_finite_int(using_gamma,
//mrb                                          nt,
//mrb                                          elecbeta,
//mrb                                          radius,
//mrb                                          Rc,
//mrb                                          expo,
//mrb                                          np,
//mrb                                          vectsfinbeta,
//mrb                                          tipo,
//mrb                                          gamma_couple,
//mrb                                          NC_minus);
//mrb        }
//mrb        
//mrb        if(radius >= Rc){
//mrb          rho = rho_radial_finite_ext(nt, 
//mrb                                      elecalfa,
//mrb                                      radius,
//mrb                                      Rc,
//mrb                                      zetas,
//mrb                                      mang,
//mrb                                      vectsfinalfa,
//mrb                                      tipo,
//mrb                                      gamma_couple,
//mrb                                      NC_plus);
//mrb          
//mrb          rho = rho + rho_radial_finite_ext(nt,
//mrb                                            elecbeta,
//mrb                                            radius,
//mrb                                            Rc,
//mrb                                            zetas,
//mrb                                            mang,
//mrb                                            vectsfinbeta,
//mrb                                            tipo,
//mrb                                            gamma_couple,
//mrb                                            NC_plus);
//mrb        }
//mrb        if (rho < rho_ini) {
//mrb          shannon = -rho*log(rho)*4.f*numero_pi*radius*radius;
//mrb          slope            = (shannon - shannon_i)/(radius - radius_i);
//mrb          intercept        = shannon_i - radius_i*slope;
//mrb          shannon_entropy = shannon_entropy + 
//mrb                            slope*(radius*radius - radius_i*radius_i)/2.f + 
//mrb                            intercept*(radius - radius_i);
//mrb          shannon_i = shannon;
//mrb          radius_i  = radius;
//mrb          rho_ini = rho;
//mrb          radius = 0.0001f + radius;
//mrb        }
//mrb        else {
//mrb          rho = 0.f;
//mrb          shannon = 0.f;
//mrb        }
//mrb       } while(fabs(rho) > bound_rho);
//mrb
//mrb       printf("\nR_inf= %6.3f, SHANNON ENTROPY UHF %f\n", radius, shannon_entropy);

       }//label 3 

  return 0;
 
 }//label 1


