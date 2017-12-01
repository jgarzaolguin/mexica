/* Integrales monoelectronicas y bielectronicas
   para resolver las ecuaciones de Hartree-Fock
   en sistemas atomicos libres.
   Jorge Garza, Junio del 2006 */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>



double bielectronic_integral_finite(char   *using_gamma,
                                    int     ang, 
				    int     mu, 
				    int     nu, 
				    int     lam, 
				    int     sig, 
				    int    *np, 
                                    int    *mang, 
				    double *zeta, 
				    double *alfa, 
				    double  Rc, 
                                    double *arreglo_factorial, 
				    double *arreglo_inv_factorial,
		                    double  gamma_couple,
                                    double *N_minus, 
				    double *N_plus)
{
 int i, a, a_1, a_2, entero_1, entero_2;
 double total, total_1, total_2, total_3, total_4, total_5, total_6,
        part_1, part_2, part_3, sum;
 double temp_1, temp_2, temp_3;
 double part_1p, part_2p, part_3p, part_1s, part_2s, part_3s;
 double temp_aux;

 double zeta_1, zeta_2, alfa_1, alfa_2;
 double zetas, alfas;
 int n_1, n_2, enes;
 int l_1, l_2, eles; 


double intc_monemln_1z_1;
double intc_mln_1z_1;
double intc_onemln_1z_1;


 extern double int1l(int arg1, double arg2, double *arreglo_factorial);
 extern double intc(int a, double b, double r, double*, double*);
 extern double upper_incomplete_gamma(double Rc, int a, double b);
 extern double use_upper_incomplete_gamma(double Rc, int a1, int a2, double b1, double b2);

 zeta_1 = zeta[mu] + zeta[nu];
 zeta_2 = zeta[lam] + zeta[sig];
 zetas = zeta_1 + zeta_2;

 alfa_1 = alfa[mu] + alfa[nu];
 alfa_2 = alfa[lam] + alfa[sig];
 alfas = alfa_1 + alfa_2;

 n_1 = np[mu] + np[nu];
 n_2 = np[lam] + np[sig];
 enes = n_1 + n_2;

 l_1 = mang[mu] + mang[nu];
 l_2 = mang[lam] + mang[sig];
 eles = l_1 + l_2;

 if (strcmp(using_gamma,"YES") == 0) {//label 1
 
    ////Constant term
     part_1 = Rc*Rc*arreglo_factorial[ang + n_2];
     part_1 = part_1/pow(zeta_2,(double) 1 + ang + n_2);
   
    ////////////////////////
    //
   
   
   
    intc_monemln_1z_1 = intc(-1 - ang + n_1, zeta_1, Rc, arreglo_factorial, arreglo_inv_factorial);
   
    intc_mln_1z_1 = intc(-ang + n_1, zeta_1, Rc, arreglo_factorial, arreglo_inv_factorial);
   
    intc_onemln_1z_1 = intc(1 - ang + n_1, zeta_1, Rc, arreglo_factorial, arreglo_inv_factorial);
   
    ///////////////////////
    //
   
    ////term1
     sum = 0.f;
     for (i = 0; i <= ang + n_2; i++)
       sum = sum +
             intc(-1 - ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
   
     sum =Rc*Rc*(intc_monemln_1z_1 - sum);
   
     temp_1 = sum;
   
    ////term2
     sum = 0.f;
     for (i = 0; i <= ang + n_2; i++)
       sum = sum + intc(-ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = -2.f*Rc*gamma_couple*(intc_mln_1z_1 - sum);
   
     temp_1 = sum + temp_1;
   
    ////term3
     sum = 0.f;
     for (i = 0; i <= ang + n_2; i++)
       sum = sum + intc(1 - ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = gamma_couple*gamma_couple*(intc_onemln_1z_1 - sum);
    
    
     temp_1 = sum + temp_1;
     temp_1 = part_1*temp_1;
    
    ////////////constant term
    
     part_1 = 0.f;
    
     part_1 = Rc*gamma_couple*arreglo_factorial[1 + ang + n_2];
     part_1 = part_1/pow(zeta_2,(double) 2 + ang + n_2);
    
    ////term4
     sum = 0.f;
     for (i = 0; i <= 1 + ang + n_2; i++)
       sum = sum + intc(-1 - ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
     sum = -2.f*Rc*Rc*(intc_monemln_1z_1 - sum);
    
     temp_2 = sum;
    
    ////term5
     sum = 0.f;
     for (i = 0; i <= 1 + ang + n_2; i++)
       sum = sum + intc(-ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
     sum = 4.f*Rc*gamma_couple*(intc_mln_1z_1 - sum);
    
     temp_2 = sum + temp_2;
    
    ///term6
     sum = 0.f;
     for (i = 0; i <= 1 + ang + n_2; i++)
       sum = sum + intc(1 - ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
     sum = -2.f*gamma_couple*gamma_couple*(intc_onemln_1z_1 - sum);
    
     temp_2 = sum + temp_2;
     temp_2 = part_1*temp_2;
    
    /////////
    /////////constant term
    
     part_1 = 0.f;
    
     part_1 = gamma_couple*gamma_couple*arreglo_factorial[2 + ang + n_2];
     part_1 = part_1/pow(zeta_2,(double) 3 + ang + n_2);
    
    ////term7
     sum = 0.f;
     for (i = 0; i <= 2 + ang + n_2; i++)
       sum = sum + intc(-1 - ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
     sum = Rc*Rc*(intc_monemln_1z_1 - sum);
    
     temp_3 = sum;
    
    ////term8
     sum = 0.f;
     for (i = 0; i <= 2 + ang + n_2; i++)
       sum = sum + intc(-ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
     sum = -2.f*Rc*gamma_couple*(intc_mln_1z_1 - sum);
    
     temp_3 = sum + temp_3;
    
    ////term9
     sum = 0.f;
     for (i = 0; i <= 2 + ang + n_2; i++)
       sum = sum + intc(1 - ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
     sum = gamma_couple*gamma_couple*(intc_onemln_1z_1 - sum);
    
     temp_3 = sum + temp_3;
     temp_3 = part_1*temp_3;
    
     part_1 = 0.f;
     part_1 = N_minus[mu]*N_minus[nu]*N_minus[lam]*N_minus[sig];
    
     total_1 = part_1*(temp_1 + temp_2 + temp_3);
     
    /////finishing first part
    //
    //jgo printf("Total 1 = %f, %f, %f, %f, %f\n", part_1, total_1, temp_1, temp_2, temp_3);
    
     part_1 = 0.f;
     part_1 = N_plus[mu]*N_plus[nu]*N_minus[lam]*N_minus[sig];
     entero_1 = l_1 + ang + 1;
     part_1 = part_1*(Rc*Rc*intc(n_2 + ang, zeta_2, Rc, arreglo_factorial, arreglo_inv_factorial) -
              2.f*Rc*gamma_couple*intc(n_2 + ang + 1, zeta_2, Rc, arreglo_factorial, arreglo_inv_factorial) +
              gamma_couple*gamma_couple*intc(n_2 + ang + 2, zeta_2, Rc, arreglo_factorial, arreglo_inv_factorial));
    /* jgo
     part_1 = part_1*upper_incomplete_gamma(Rc, entero_1, alfa_1);
    jgo */
     part_1 = 0.f;
     total_2 = part_1;
    
    /////////////
    
     entero_2 = mang[lam] + mang[sig] - ang;
    /* jgo
     part_1 = upper_incomplete_gamma(Rc, entero_1, alfa_1);
     part_1 = part_1*upper_incomplete_gamma(Rc, entero_2, alfa_2);
     part_1 = part_1 - use_upper_incomplete_gamma(Rc, entero_1,  entero_2, alfa_1, alfa_2);
     part_1 = part_1*N_plus[mu]*N_plus[nu]*N_plus[lam]*N_plus[sig];
    jgo */
     part_1 = 0.f;
     total_3 = part_1;
    
    //////////
    
     temp_aux = arreglo_factorial[-1 - ang + n_2];
     temp_aux = temp_aux/pow(zeta_2,(double) -ang + n_2);
    
    ////term1
    
     part_1p = intc(ang + n_1, zeta_1, Rc, arreglo_factorial, arreglo_inv_factorial);
     part_1s = intc(-1 - ang + n_2, zeta_2, Rc, arreglo_factorial, arreglo_inv_factorial);
    
     part_1 = 0.f;
     part_1 = part_1p*part_1s;
    
    
     sum = 0.f;
     for (i = 0; i <= -1 - ang + n_2; i++)
       sum = sum + intc(ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = Rc*Rc*Rc*Rc*(temp_aux*(sum - part_1p) + part_1);
    
     temp_1 = sum;
    
    
    
    ////term2
    
     part_2p = intc(1 + ang + n_1, zeta_1, Rc, arreglo_factorial, arreglo_inv_factorial);
     
     part_1 = 0.f;
     part_1 = part_2p*part_1s;
    
      sum = 0.f;
     for (i = 0; i <= -1 - ang + n_2; i++)
       sum = sum + intc(1 + ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = -2.f*Rc*gamma_couple*Rc*Rc*(temp_aux*(sum - part_2p) + part_1);
    
     temp_1 = temp_1 + sum;
    
    ////term3
    
     part_3p = intc(2 + ang + n_1, zeta_1, Rc, arreglo_factorial, arreglo_inv_factorial);
    
      part_1 = 0.f;
      part_1 = part_3p*part_1s;
    
      sum = 0.f;
     for (i = 0; i <= -1 - ang + n_2; i++)
       sum = sum + intc(2 + ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = gamma_couple*gamma_couple*Rc*Rc*(temp_aux*(sum - part_3p) + part_1);
    
     
     temp_1 = temp_1 + sum;
    
    ////term4
    
     temp_aux = 0.f;
     temp_aux = arreglo_factorial[-ang + n_2];
     temp_aux = temp_aux/pow(zeta_2,(double) 1 - ang + n_2);
    
      
     part_2s = intc(-ang + n_2, zeta_2, Rc, arreglo_factorial, arreglo_inv_factorial);
    
     part_1 = 0.f;
     part_1 = part_1p*part_2s;
    
       sum = 0.f;
     for (i = 0; i <= -ang + n_2 ; i++)
       sum = sum + intc(ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = -2.f*Rc*Rc*Rc*gamma_couple*(temp_aux*(sum - part_1p) + part_1);
    
     temp_2 = sum;
    
    
    ////term5
    
     part_1 = 0.f;
     part_1 = part_2p*part_2s;
    
       sum = 0.f;
     for (i = 0; i <= -ang + n_2 ; i++)
       sum = sum + intc(1 + ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = 4.f*Rc*gamma_couple*Rc*gamma_couple*(temp_aux*(sum - part_2p) + part_1);
    
     temp_2 = temp_2 + sum;
     
    ////term6
    
     part_1 = 0.f;
     part_1 = part_3p*part_2s;
    
       sum = 0.f;
     for (i = 0; i <= -ang + n_2 ; i++)
       sum = sum + intc(2 + ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = -2.f*gamma_couple*gamma_couple*Rc*gamma_couple*(temp_aux*(sum - part_3p) + part_1);
    
     temp_2 = temp_2 + sum;
    
    ////term7
    
     temp_aux = 0.f;
     temp_aux = arreglo_factorial[1 - ang + n_2];
     temp_aux = temp_aux/pow(zeta_2,(double) 2 - ang + n_2);
    
     part_3s = intc(1 - ang + n_2, zeta_2, Rc, arreglo_factorial, arreglo_inv_factorial);
    
     part_1 = 0.f;
     part_1 = part_1p*part_3s;
    
        sum = 0.f;
     for (i = 0; i <= 1 - ang + n_2; i++)
       sum = sum + intc(ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = Rc*Rc*gamma_couple*gamma_couple*(temp_aux*(sum - part_1p) + part_1);
    
     temp_3 = sum;
    
    ////term8
    
      part_1 = 0.f;
      part_1 = part_2p*part_3s;
    
        sum = 0.f;
     for (i = 0; i <= 1 - ang + n_2; i++)
       sum = sum + intc(1 + ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = -2.f*Rc*gamma_couple*gamma_couple*gamma_couple*(temp_aux*(sum - part_2p) + part_1);
    
     temp_3 = temp_3 + sum;
    
    ////term9
    
      part_1 = 0.f;
      part_1 = part_3p*part_3s;
    
        sum = 0.f;
     for (i = 0; i <= 1 - ang + n_2; i++)
       sum = sum + intc(2 + ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
             pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
    
     sum = gamma_couple*gamma_couple*gamma_couple*gamma_couple*(temp_aux*(sum - part_3p) + part_1);
    
     temp_3 = temp_3 + sum;
    
    
     temp_aux = 0.f;
     temp_aux = temp_1 + temp_2 + temp_3;
    
     temp_aux = N_minus[mu]*N_minus[nu]*N_minus[lam]*N_minus[sig]*temp_aux;
    
     total_4 = temp_aux;
    
    //////////////////////////////
    //////////////////////////////
    /////////////////////////////
    
     part_1 = 0.f;
     part_1 = N_minus[mu]*N_minus[nu]*N_plus[lam]*N_plus[sig];
     part_1 = part_1*(Rc*Rc*part_1p - 2.f*Rc*gamma_couple*part_2p + gamma_couple*gamma_couple*part_3p);
    /* jgo
     part_1 =part_1*upper_incomplete_gamma(Rc, l_2 + ang + 1, alfa_2);
    jgo */
     part_1 = 0.f;
    
     total_5 = part_1;
    
     part_1 = 0.f;
    // part_1 = use_upper_incomplete_gamma(Rc, mang[mu] + mang[nu] - ang,  mang[lam] + mang[sig] + ang + 1, alfa_1, alfa_2);
    /* jgo
     part_1 = use_upper_incomplete_gamma(Rc, l_1 - ang,  l_2 + ang + 1, alfa_1, alfa_2);
     part_1 = part_1*N_plus[mu]*N_plus[nu]*N_plus[lam]*N_plus[sig];
    jgo */
     part_1 = 0.f;
    
     total_6 = part_1;
    // printf("jgo 4: totales = %f, %f, %f, %f, %f, %f\n", total_1, total_2, total_3, total_4, total_5, total_6);
     total = total_1 + total_2 + total_3 + total_4 + total_5 + total_6;
   }//label 1
    else {//label 2

         a = ang + n_2;
         
         part_1 = arreglo_factorial[a];
         part_1 = part_1/pow(zeta_2,(double) a + 1);
         
         sum = 0.f;
         for (i = 0; i <= a; i++)
           sum = sum +
                 intc(-1 - ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*
                 pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
         
         sum = intc(-1 - ang + n_1, zeta_1, Rc, arreglo_factorial, arreglo_inv_factorial) - sum;
         
         total_1 = part_1*N_minus[mu]*N_minus[nu]*N_minus[lam]*N_minus[sig]*sum;
         
         
         part_1 = 0.f;
         /* mrb
          part_1 = N_plus[mu]*N_plus[nu]*N_minus[lam]*N_minus[sig];
          entero_1 = l_1 + ang + 1;
          part_1 = part_1*intc(n_2 + ang, zeta_2, Rc, arreglo_factorial, arreglo_inv_factorial);
          part_1 = part_1*upper_incomplete_gamma(Rc, entero_1, alfa_1);
         */
         total_2 = part_1;
         
         entero_2 = l_2 - ang;
         
         part_1 = 0.f;
         /*mrb
          part_1 = upper_incomplete_gamma(Rc, entero_1, alfa_1);
          part_1 = part_1*upper_incomplete_gamma(Rc, entero_2, alfa_2);
          part_1 = part_1 - use_upper_incomplete_gamma(Rc, entero_1,  entero_2, alfa_1, alfa_2);
          part_1 = part_1*N_plus[mu]*N_plus[nu]*N_plus[lam]*N_plus[sig];
         */
         total_3 = part_1;
         
         a = -1 - ang + n_2;
         
         temp_aux = arreglo_factorial[a];
         temp_aux = temp_aux/pow(zeta_2,(double) a + 1);
         
         part_1p = intc(ang + n_1, zeta_1, Rc, arreglo_factorial, arreglo_inv_factorial);
         part_1s = intc(-1 - ang + n_2, zeta_2, Rc, arreglo_factorial, arreglo_inv_factorial);
         
         part_1 = 0.f;
         part_1 = part_1p*part_1s;
         
         
         sum = 0.f;
         for (i = 0; i <= a; i++)
           sum = sum + intc(ang + i + n_1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial)*pow(zeta_2,(double) i)*arreglo_inv_factorial[i];
         
         sum = temp_aux*(sum - part_1p) + part_1;
         
         total_4 = sum*N_minus[mu]*N_minus[nu]*N_minus[lam]*N_minus[sig];
         
         
         part_1 = 0.f;
         /* mrb
          part_1 = N_minus[mu]*N_minus[nu]*N_plus[lam]*N_plus[sig];
          part_1 = part_1*intc(n_1 + ang, zeta_1, Rc, arreglo_factorial, arreglo_inv_factorial);
          part_1 = part_1*upper_incomplete_gamma(Rc, l_2 + ang + 1, alfa_2);
         */
         
         total_5 = part_1;
         
         part_1 = 0.f;
         /* mrb
          part_1 = use_upper_incomplete_gamma(Rc, l_1 - ang,  l_2 + ang + 1, alfa_1, alfa_2);
          part_1 = part_1*N_plus[mu]*N_plus[nu]*N_plus[lam]*N_plus[sig];
         */
         
         total_6 = part_1;
         total = total_1 + total_2 + total_3 + total_4 + total_5 + total_6;


        }//label 2


   return (total);
  }
  
