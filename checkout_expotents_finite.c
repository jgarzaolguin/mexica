#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


	int check_expotents_finite(char   *using_gamma,
                                   double  gamma_couple,
				   double  Rc,
				   int    *np,
			           int    *mang,
			 	   double *expo,
		  		   int     index)
	{
          int    result;
	  double factor1,
	         factor2,
		 alpha_mu,
                 temp_2,
                 exponente_interno_temporal;
    
       if (strcmp(using_gamma,"YES") == 0){//label 1

          factor1  = gamma_couple/((1.f - gamma_couple)*Rc) + expo[index];
          factor2  = (double) (mang[index] + np[index])/Rc;
          alpha_mu = factor1 - factor2;

          result = 0;
          if (factor1 < factor2 || alpha_mu <= 0.01 || expo[index] <= 0.01)
             result = 1;
          }//label 1
           else {//label 2
          

           temp_2   = (double) (np[index] + mang[index]);
           exponente_interno_temporal = expo[index];
           alpha_mu = exponente_interno_temporal - temp_2/Rc;

           result = 0;
           if (exponente_interno_temporal < temp_2/Rc || alpha_mu <= 0.01 || exponente_interno_temporal <= 0.01)
             result = 1;

               }//label 2




          return result;
	}



