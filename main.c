#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
// MEXICA-C code to study free atoms or atoms confined
// by penetrable or impenetrable walls.
// Authors: Jorge Garza (jgo@xanum.uam.mx)
//          Mariano Rodriguez-Bautista (marianorodbau@gmail.com)
// Version 1.0. Sept, 2015.
// Please cite this work as:
//--------------------------------------------------------------------------------------------------
//            Jorge Garza, Julio M. Hern\'ande-P\'erez, Jos\'e-Zeferino Ram\'irez and Rubicelia Vargas.
//            Basis set effects on the Hartree-Fock description of confined many-electron atoms.
//            Journal of Physics B: Atomic, Molecular & Optical Physics. 45, 015002 (2012)
//--------------------------------------------------------------------------------------------------
//            Mariano Rodriguez-Bautista, Cecilia D\'iaz-Garc\'ia, Alejandra M. Navarrete-L\'opez,
//            Rubicelia Vargas, and Jorge Garza.
//            Roothaan's approach to solve the Hartree-Fock equations for atoms 
//            confined by soft walls: Basis set with correct asymptotic behavior.
//            The Journal of Chemical Physics. 143, 034103 (2015).
//--------------------------------------------------------------------------------------------------

int optimiza_main(int     nt, 
                  char   *opt,
		  int     elecalfa, 
		  int     elecbeta,
                  double  z, 
		  int     orbital, 
		  double  tol, 
                  char   *using_gamma,
                  char   *correlation,
                  char   *propagador,
                  char   **save_dft,
                  double *weight_dft,
                  int     flag_dft,
		  double  mezcla, 
                  char   *basis,
		  double *expo, 
		  int    *np,
                  int    *mang, 
		  int    *ncm,  
		  double  gamma_couple,
                  char   *tipo, 
		  int     maxiter, 
		  double  Rc, 
		  char   *bound,
                  char   *espin, 
		  double  step, 
		  int    *opt_flag,
                  double  epsilon, 
		  char   *config,
                  char   *nombre, 
		  int     imprime,
                  int plasma,
                  char   *properties)
{ //Function to optimize exponents in the basis set
 int i, j;
 double maxdiff;
 FILE* archivo_opt;

// double expo_diff[800];

 extern void optimiza(int     nt, 
                      char   *opt,
		      int     elecalfa, 
		      int     elecbeta,
                      double  z, 
		      int     orbital, 
		      double  tol, 
                      char   *using_gamma,
                      char   *correlation,
                      char   *propagador,
                      char   **save_dft,
                      double *weight_dft,
                      int     flag_dft,
		      double  mezcla, 
                      char   *basis,
		      double *expo, 
		      int    *np,
                      int    *mang, 
		      int    *ncm, 
		      double  gamma_couple,
                      char   *tipo, 
		      int     maxiter, 
		      double  Rc, 
                      char   *nombre,
		      char   *bound,
                      char   *espin, 
		      double  step, 
		      int    *opt_flag,
                      double  epsilon, 
		      int     imprime,
                      int plasma,
                      char   *properties);

//mrb    for (i = 0; i < nt; i++) expo_diff[i] = 0.f;
//mrb    i = 0;

//mrb    do {
//mrb         maxdiff = 0.f;
//mrb         i = i + 1;
//mrb         printf("Step %d in the optimization process\n", i);

          optimiza(nt, 
                   opt,
		   elecalfa, 
	    	   elecbeta, 
		   z, 
		   orbital, 
		   tol, 
                   using_gamma,
                   correlation,
                   propagador,
                   save_dft,
                   weight_dft,
                   flag_dft,
		   mezcla,
                   basis,
                   expo, 
		   np, 
		   mang, 
		   ncm, 
	    	   gamma_couple, 
		   tipo, 
                   maxiter, 
		   Rc, 
                   nombre,
		   bound, 
		   espin, 
		   step, 
	  	   opt_flag,
                   epsilon, 
		   imprime,
                   plasma,
                   properties);

//mrb          for (j = 0; j < nt; j++) 
//mrb             {
//mrb                 if(fabs(expo[j] - expo_diff[j]) > maxdiff)
//mrb                    maxdiff = fabs(expo[j] - expo_diff[j]);
//mrb                 expo_diff[j] = expo[j];
//mrb             }
//mrb       } while (maxdiff > 1.e-5);

     printf("Reading initial basis set from %s\n", nombre);
     archivo_opt = fopen(nombre,"r");
     printf(":)  Final exponents:\n");
     double temp_dble;
     int temp_int;
     i = -1;
 printf("\nENTRAAAA\n");
    
     while (fscanf(archivo_opt,"%s %lf %d", config, &temp_dble, &temp_int) != EOF)
     {
       if (config[1] == 'S') i = i + 1;
       if (config[1] == 'P') i = i + 3;
       if (config[1] == 'D') i = i + 5;
       if (config[1] == 'F') i = i + 7;
       if (config[1] == 'G') i = i + 9;
       if (config[1] == 'H') i = i + 11;
       if (config[1] == 'I') i = i + 13;
       printf(":)  %s  %8.25lf %d \n", config, expo[i], temp_int);
     }
     fclose(archivo_opt);
     printf("** SCF after optimization **\n");
 return (0);
 }

 int main(int argc, char *argv[])
 {//Main function, this functions controls SCF and optimization of exponents
  //in the basis set
  int i, j, test_inp;
  int nt, elecalfa, elecbeta, maxiter, orbital;
  int prue1, prue2, prue3, prue4, plasma;
 
  double energia, z, tol, Rc, step, mezcla, maxdiff, epsilon, gamma_couple;
  double x0, x1, x2, energia_0, energia_1, temp_gam;
  char config[4], read_base[10], bound[15], tipo[16],  opt[10],
       espin[10], nombre[100], correlation[10], propagador[10], using_gamma[10],
       properties[10], basis[15];
  double energy_array[500], gamma_array[500];
 
  char kind_of_cal[10];
  double count_temp, count_final, add_tol;
  int steps;

  double exponente_interno_temporal,
         factor1,
         factor2;
 
  char   **save_dft;

      save_dft = (char **) malloc(80*sizeof(char *));
   if (save_dft == NULL) {
     free(save_dft);
     save_dft = 0;
     printf("Error en arreglo save_dft, insuficiente RAM\n");
     exit(1);
   }

for(i=0 ; i < 80 ; i++){
    save_dft[i] = (char *) malloc(10*sizeof(char));
    if (save_dft[i] == NULL) {
     free(save_dft[i]);
     save_dft[i] = 0;
     printf("Error en arreglo save_dft, insuficiente RAM\n");
     exit(1);
   }
}



  double weight_dft[80];
 for (i = 0; i < 80; i++){
    strcpy(save_dft[i],"end");
    weight_dft[i] = 0.f;
 }


 int  cont;
  
 FILE* archivo_opt;
 
 extern int         input(double   *z, 
 	                    int    *elecalfa, 
 	                    int    *elecbeta, 
                            char   *tipo, 
                            char   *using_gamma,
                            char   *correlation,
                            char   *propagador,
                            char   **save_dft,
                            double *weight_dft,
                            int    *flag_dft,
                            char   *bound, 
                            char   *espin,
                            double *Rc, 
                            double *tol, 
                            double *mezcla, 
                            int    *orbital,
                            int    *maxiter, 
                            char   *basis,
                            char   *read_base,
                            int    *nt, 
                            int    *np, 
                            int    *mang, 
                            int    *ncm, 
                            double *gamma_couple,
                            double *expo, 
                            char   *opt, 
                            double *step, 
                            int    *opt_flag,
                            char   *nombre, 
                            double *epsilon, 
                            char   *kind_of_cal, 
                            double *count_temp, 
                            double *count_final, 
                            int    *steps,
                            int    *plasma,
                            char   *properties);
  
 extern int           scf(int     nt, 
                           int     elecalfa, 
                           int     elecbeta,
                           double  z, 
                           int     orbital,
                           double  tol, 
                           char   *using_gamma,
                           char   *correlation,
                           char   *propagador,
                           char   **save_dft,
                           double  *weight_dft,
                           int      flag_dft,
                           double  mezcla,
                           char   *basis,
                           double *expo, 
                           int    *np, 
                           int    *mang, 
                           int    *ncm, 
                           double  gamma_couple,
                           char   *tipo, 
                           int     maxiter, 
                           double  Rc,
                           char   *bound, 
                           char   *espin, 
                           double *total_energy, 
                           int     print_vectors, 
                           double  epsilon,
                           int     imprime,
                           int     plasma,
                           double *cusp_kato,
                           char   *properties);
 
 extern int optimiza_main( int     nt, 
                           char   *opt,
                           int     elecalfa, 
                           int     elecbeta,
                           double  z, 
                           int     orbital, 
                           double  tol, 
                           char   *using_gamma,
                           char   *correlation,
                           char   *propagador,
                           char   **save_dft,
                           double  *weight_dft,
                           int      flag_dft,
                           double  mezcla, 
                           char   *basis,
                           double *expo, 
                           int    *np,
                           int    *mang, 
                           int    *ncm, 
                           double  gamma_couple,
                           char   *tipo, 
                           int     maxiter, 
                           double  Rc, 
                           char   *bound,
                           char   *espin, 
                           double  step, 
                           int    *opt_flag,
                           double  epsilon, 
                           char   *config,
                           char   *nombre, 
                           int     imprime,
                           int     plasma,
                           char   *properties);

 extern int  main_finite_global(double   z,
                               int      elecalfa,
                               int      elecbeta,
                               char    *tipo,
                               char    *using_gamma,
                               char    *correlation,
                               char    *propagador,
                               char   **save_dft,
                               double  *weight_dft,
                               int      flag_dft,
                               char    *bound,
                               char    *espin,
                               double   Rc,
                               double   tol,
                               double   mezcla,
                               int      orbital,
                               int      maxiter,
                               int      nt,
                               int     *np,
                               int     *mang,
                               int     *ncm,
                               double   gamma_couple,
                               double  *expo,
                               char    *opt,
                               double   step,
                               int     *opt_flag,
                               char    *nombre,
                               double   epsilon,
                               char    *kind_of_cal,
                               double   count_temp,
                               double   count_final,
                               int      steps,
                               double  *energia,
                               int      plasma,
                               char    *properties);


 
 // Maximum size of the Fock matrix and related matrices, you can change it and recompile
 int maxfun = 800;
 int np[maxfun], mang[maxfun], ncm[maxfun], opt_flag[maxfun];
 double expo[maxfun];
 double cusp_kato;
 int    flag_dft;  
 // Initializing arrays associated to the basis set
 for (i = 0; i < maxfun; i++) {
    np[i] = 0;
    mang[i] = 0;
    ncm[i] = 0;
    expo[i] = (double) 0.;
 }
 
 test_inp = input(&z, 
 		  &elecalfa, 
 	   	  &elecbeta, 
 		   tipo, 
                   using_gamma,
                   correlation,
                   propagador,
                   save_dft,
                   weight_dft,
                  &flag_dft,
                   bound, 
  		   espin, 
 		  &Rc, 
 	  	  &tol, 
  		  &mezcla, 
 		  &orbital, 
 		  &maxiter, 
                   basis,
                   read_base, 
 		  &nt, 
 		   np, 
 	  	   mang, 
  		   ncm, 
 		  &gamma_couple, 
 		   expo, 
 		   opt, 
 		  &step,
                   opt_flag, 
  		   nombre, 
 		  &epsilon,
                   kind_of_cal, 
 		  &count_temp, 
  		  &count_final, 
 		  &steps,
                  &plasma,
                   properties);

// for (i = 0; i < 80; i++) {
// printf("%s %lf %d\n", save_dft[i], weight_dft[i], flag_dft);
// }

  
 if(test_inp == 0) {
    if(strcmp(opt,"opt") == 0 || strcmp(opt,"optfull") == 0) {       //Optimization of exponents in the basis set
       if(strcmp(basis,"STOs") == 0 || strcmp(basis,"stos") == 0) {  //begins stos
          if(strcmp(bound,"finite") == 0 || strcmp(bound,"dielectricc") == 0 || strcmp(bound,"polarization") == 0) { //Confinement imposed by penetrable walls
             main_finite_global(z,
                                elecalfa,
                                elecbeta,
                                tipo,
                                using_gamma,
                                correlation,
                                propagador,
                                save_dft,
                                weight_dft,
                                flag_dft,
                                bound,
                                espin,
                                Rc,
                                tol,
                                mezcla,
                                orbital,
                                maxiter,
                                nt,
                                np,
                                mang,
                                ncm,
                                gamma_couple,
                                expo,
                                opt,
                                step,
                                opt_flag,
                                nombre,
                                epsilon,
                                kind_of_cal,
                                count_temp,
                                count_final,
                                steps,
                                &energia,
                                plasma,
                                properties);
          } 
          else { //Procedure for free atoms or confinement imposed by impenetrable walls
             optimiza_main(nt, 
                           opt,
                           elecalfa, 
                           elecbeta,
                           z,  
                           orbital,  
                           tol,  
                           using_gamma,
                           correlation,
                           propagador,
                           save_dft,
                           weight_dft,
                           flag_dft,
                           mezcla, 
                           basis,
                           expo, 
                           np,
                           mang, 
                           ncm, 
                           gamma_couple,
                           tipo, 
                           maxiter, 
                           Rc, 
                           bound,
                           espin, 
                           step, 
                           opt_flag,
                           epsilon, 
                           config, 
                           nombre, 
                           0,
                           plasma,
                           properties);
             scf(nt, 
                 elecalfa, 
                 elecbeta, 
                 z, 
                 orbital, 
                 tol, 
                 using_gamma,
                 correlation,
                 propagador,
                 save_dft,
                 weight_dft,
                 flag_dft,
                 mezcla, 
                 basis,
                 expo, 
                 np, 
                 mang,
                 ncm, 
                 gamma_couple, 
                 tipo, 
                 maxiter, 
                 Rc, 
                 bound, 
                 espin,
                 &energia, 
                 1, 
                 epsilon, 
                 1,
                 plasma,
                 &cusp_kato,
                 properties);
          }
       }    // ends stos
       else {
          if(strcmp(basis,"GTOs") == 0 || strcmp(basis,"gtos") == 0) {
             optimiza_main(nt,
                           opt,
                           elecalfa,
                           elecbeta,
                           z,
                           orbital,
                           tol,
                           using_gamma,
                           correlation,
                           propagador,
                           save_dft,
                           weight_dft,
                           flag_dft,
                           mezcla,
                           basis,
                           expo,
                           np,
                           mang,
                           ncm,
                           gamma_couple,
                           tipo,
                           maxiter,
                           Rc,
                           bound,
                           espin,
                           step,
                           opt_flag,
                           epsilon,
                           config,
                           nombre,
                           0,
                           plasma,
                           properties); 
             scf(nt,
                 elecalfa,
                 elecbeta,
                 z,
                 orbital,
                 tol,
                 using_gamma,
                 correlation,
                 propagador,
                 save_dft,
                 weight_dft,
                 flag_dft,
                 mezcla,
                 basis,
                 expo,
                 np,
                 mang,
                 ncm,
                 gamma_couple,
                 tipo,
                 maxiter,
                 Rc,
                 bound,
                 espin,
                 &energia,
                 1,
                 epsilon,
                 1,
                 plasma,
                 &cusp_kato,
                 properties);
          }
       }    // ends gtos
    }       // ends opt calculation 
    else {  //Just one SCF calculation
       scf(nt, 
           elecalfa, 
           elecbeta, 
           z, 
           orbital, 
           tol, 
           using_gamma,
           correlation,
           propagador,
           save_dft,
           weight_dft,
           flag_dft,
           mezcla, 
           basis,
           expo, 
           np, 
           mang,
           ncm, 
           gamma_couple, 
           tipo, 
           maxiter, 
           Rc, 
           bound, 
           espin,
           &energia, 
           1, 
           epsilon, 
           1,
           plasma,
           &cusp_kato,
           properties);
    }
//
// Reporting final energy
//
    printf("------------------------------------\n");
    printf("Calculation done\n");
    printf("Energy %15.5f a.u  \n",energia);
    printf("Energy %15.5f Ryd\n", 2.f*energia);
    printf("Energy %15.5f eV \n", 27.211399*energia);
//    printf("Energy %15.5f kJ/mol \n", 2625.5002*energia);
//    printf("Energy %15.5f Kcal/mol\n",627.5096080305927*energia);
//printf("Energy %15.5f K (equivalent temperature)\n",315777.09*energia);
//    printf("Energy %15.5f cm**-1\n", 219474.63*energia);
//ADRIAN AGREGO    printf("Energy %15.5f GHz \n", 6579683.879634054*energia);
//ADRIAN AGREGO    printf("Energy %15.5f nm \n", 0.021947463*energia);
//ADRIAN AGREGO    printf("------------------------------------\n");
//ADRIAN AGREGO 
    remove(nombre);
 
   }

   for (i = 79; i <= 0; i--) free(save_dft[i]);
   free(save_dft);


  return 0;
 }
