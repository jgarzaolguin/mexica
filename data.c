// Input data management
//
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>

int input(double *z, 
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
          char   *read_base,
          int    *base, 
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
          char   *properties)
{
 FILE* file_1;
 int i, j, c, d, entero1, entero2, z_temp, charge, mult, alfa, beta, prev, test_input;
 double Rc_temp, temp, epsilon_temp, temp1, temp2, temp_g;
 int temp3;

 int sizeint, sizedouble;
 char elemento[3], oneword[80], twoword[80], threeword[3];
 char scratch[100];

 time_t t;
 struct tm tstruct;

 extern	int input_2_int(char   *read_base,
       		        int     base,
                        double  Rc,
       		        int    *np_temp,
       		        int    *mang_temp,
       		        int    *ncm_temp,
       		        double *expo_temp,
       		        int    *nt,
       		        char   *opt,
       		        int    *opt_flag,
       		        FILE   *file_1,
       		        char   *nombre);

 extern	int input_2_ext(char   *read_base,
       		        int     base,
                        double  Rc,
       		        int    *np_temp,
       		        int    *mang_temp,
       		        int    *ncm_temp,
       		        double *expo_temp,
       		        int    *nt,
       		        int    *opt_flag,
       		        FILE   *file_1,
       		        char   *nombre);

 int atomic_number[54] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                  11, 12, 13, 14, 15, 16, 17, 18,
                  19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                  31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
                  43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 52, 54};

 char symbol[54][3] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
           "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
           "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
           "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
           "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
           "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"};

 test_input = 0;

 printf("----------------------------------\n");
 printf("-- MEXICA \n");

 t = time(0);
 tstruct = *localtime(&t);
 printf("-- %s \n", ctime(&t));

 scanf("%d", &z_temp);
 *z = (double)z_temp;     // Nuclear charge

 scanf("%d %d", &charge, &mult);   // Charge and multiplicity
 prev = mult + z_temp - charge - 1; 
 alfa = prev/2;
 beta = z_temp - charge - alfa;
   
 if ( fmod(prev,2) != (double)0. || beta < 0) {
     printf("-- ERROR -- Check charge and multiplicity --\n");
     return 1;
 } else {
     *elecalfa = alfa;
     *elecbeta = beta;
     scanf("%s %s %s", tipo, correlation, propagador); // Method to be used, with or without correlation or electron propagator
     if (strcmp(tipo,"rhf") == 0 || 
	 strcmp(tipo,"uhf") == 0 || 
	 strcmp(tipo,"rks") == 0 || 
	 strcmp(tipo,"uks") == 0 ) {
       if (mult != 1) {
         if ( strcmp(tipo,"rhf") == 0 ){
           printf("Multiplicity = %d, changing method to uhf\n", mult);
           strcpy(tipo,"uhf");
         } else {
             if ( strcmp(tipo,"rks") == 0 ){
               printf("Multiplicity = %d, changing method to uks\n", mult);
               strcpy(tipo,"uks");
             } 
         }
       }
       if (strcmp(correlation,"mp2") == 0 || strcmp(correlation,"MP2") == 0)
         strcpy(correlation,"mp2");
       if (strcmp(propagador,"ep2") == 0 || strcmp(propagador,"EP2") == 0) {
         strcpy(propagador,"tom");
         printf("Electron propagator to second order\n");
         printf("Provide alpha or beta component and orbital to detach one electron (start from zero): \n");
         scanf("%s %d", espin, orbital);
       }
     }
     else {
       printf("Check the calculation, no %s method\n", tipo);
       return 1;
     }


      //propuesta para leer funcionales dft
      scanf("%s", save_dft[0]);
      weight_dft[0] = 0.f;
      int unsigned dft; 
      int unsigned ident;
      ident = 0;
      dft = 1;

      if (strcmp(save_dft[0], "xc") == 0) {

        do {
          scanf("%s", save_dft[dft]);

          if (strcmp(save_dft[dft],"end") == 0){
             ident = 1;
             weight_dft[dft] = 0.f;
           } else {
              scanf("%lf",&weight_dft[dft]);
             }
          //printf("\nHOLA    %s %lf\n", save_dft[dft], weight_dft[dft]);
          dft++; 
        } while (ident != 1); 

     }

     *flag_dft = (dft - 1);

      scanf("%lf", &temp);
     *tol = temp;              // SCF tolerance
     
     scanf("%lf", &temp);
     *mezcla = temp;          // Mix of density matrix in the SCF

     scanf("%d", &entero1);
     *maxiter = entero1;      // Iterations for the SCF

     scanf("%s %lf", opt, &temp);
     *step = temp;            // Exponents optimization starting the changes in step
     strcpy(elemento, symbol[z_temp - 1]);
     printf("Atom : %s\n", elemento);
     printf("Alpha and beta : %d, %d\n", alfa, beta);
     Rc_temp = 1e20;
     epsilon_temp = 0e00;
     *plasma = 0;
     scanf("%s", bound);
     if (strcmp(bound,"free") == 0 || strcmp(bound,"Free") == 0 || strcmp(bound,"FREE") == 0) {
       *Rc = 0.00000;
       printf("Free atom calculation\n");
     } else
     if (strcmp(bound,"dielectricnc") == 0 || strcmp(bound,"Dielectricnc") == 0 || strcmp(bound,"DIELECTRICNC") == 0) {
       scanf("%lf %lf", &Rc_temp, &epsilon_temp);
       *Rc = Rc_temp;
       *epsilon = epsilon_temp;
       printf("Dielectric confined with asymtotic behavior incorrect\n");
     } else 
     if (strcmp(bound,"dielectricc")  == 0  || 
         strcmp(bound,"Dielectricc")  == 0  || 
         strcmp(bound,"DIELECTRICC")  == 0  || 
         strcmp(bound,"polarization") == 0  || 
         strcmp(bound,"Polarization") == 0  || 
         strcmp(bound,"POLARIZATION") == 0) {
       scanf("%lf %lf", &Rc_temp, &epsilon_temp);
       *Rc = Rc_temp;
       *epsilon = epsilon_temp;
       scanf("%s",using_gamma);
       if (strcmp(using_gamma,"n") == 0 || strcmp(using_gamma,"no") == 0 || strcmp(using_gamma,"NO") == 0 || strcmp(using_gamma,"No") == 0)
       strcpy(using_gamma,"NO");
       if (strcmp(bound,"polarization") == 0)
       printf("Polarization with asymtotic behavior correct\n");
       else
       printf("Dielectric confined with asymtotic behavior correct\n");
     }  
     else {
     if (strcmp(bound,"impenetrable") == 0 || strcmp(bound,"Impenetrable") == 0 || strcmp(bound,"IMPENETRABLE") == 0) {
       strcpy(bound,"confined");
       printf("Atom confined by impenetrable walls, please give the confinement radius\n");
       scanf("%lf", &Rc_temp);
       *Rc = Rc_temp;
       printf("Impenetrable walls at Rc=%f\n", *Rc);
     } else {
     if (strcmp(bound,"penetrable") == 0 || strcmp(bound,"Penetrable") == 0 || strcmp(bound,"PENETRABLE") == 0) {
       strcpy(bound,"finite");
       printf("Atom confined by penetrable walls,\n");
       printf("provide the confinement radius and the height of the barrier\n");
       scanf("%lf %lf", &Rc_temp, &epsilon_temp);
       *Rc = Rc_temp;
       *epsilon = epsilon_temp;
       printf("Penetrable walls at Rc=%f and height=%f\n", *Rc, *epsilon);
//     Default values
       strcpy(using_gamma,"NO");
       strcpy(kind_of_cal,"NO");
       temp_g = 7e-1;
       *gamma_couple = temp_g;
       temp1 = 1e-9;
       temp2 = 98e-2;
       temp3 = 20;
       *count_temp = temp1; 
       *count_final = temp2; 
       *steps = temp3;
       printf("Auxiliar function in the basis set? (y/n)\n");
       scanf("%s", using_gamma);
       if (strcmp(using_gamma,"y") == 0 || strcmp(using_gamma,"yes") == 0 || strcmp(using_gamma,"YES") == 0) {
         strcpy(using_gamma,"YES");
         printf("Optimization of gamma in auxiliar function? (y/n)\n");
         scanf("%s", kind_of_cal);
         if (strcmp(kind_of_cal,"y") == 0 || strcmp(kind_of_cal,"yes") == 0 || strcmp(kind_of_cal,"YES") == 0) {
           strcpy(kind_of_cal,"YES");
           temp_g = 5e-01;
           *gamma_couple = temp_g;
           printf("Gamma optimization\n");
         } else {
           printf("No optimization for gamma, please provide its value\n");
           scanf("%lf", &temp_g);
           *gamma_couple = temp_g;
           printf("Auxiliar functions without optimization of gamma, gamma=%f\n", *gamma_couple);
         }
       }
       else {
         strcpy(using_gamma,"NO");
         printf("No auxiliar funtion in the basis set\n");
       }
     } else {
         if (strcmp(bound,"plasma") == 0 || strcmp(bound,"Plasma") == 0 || strcmp(bound,"PLASMA") == 0) {
           strcpy(bound,"confined");
           printf("Atom confined by plasma model, please give the confinement radius\n");
           scanf("%lf", &Rc_temp);
           *Rc = Rc_temp;
           *plasma = 1;
           printf("Plasma model with impenetrable walls at Rc=%f\n", *Rc);
         }
     }
     }
     }

     scanf("%s", properties);
     printf("Exponents optimization : ");
     if ((strcmp(opt,"opt") == 0) || (strcmp(opt,"OPT") == 0) || (strcmp(opt,"Opt") == 0))
       printf("yes\n");
     else
       printf("no\n");

     scanf("%s", read_base);

     if (strcmp(read_base,"bunge")    == 0 || 
	 strcmp(read_base,"clementi") == 0 || 
	 strcmp(read_base,"thakkar")  == 0){
       if ( z_temp < 1 || z_temp > 54) {
         printf("Atomic Number out of range in the internal basis set.\n");
         printf("You must provide the basis set for this atomic number.\n");
         return 1;
       } else {
           if (strcmp(read_base,"bunge") == 0) 
             file_1 = fopen("base_bunge.txt", "r");
           else
           if (strcmp(read_base,"clementi") == 0)
             file_1 = fopen("base_clementi.txt", "r");
           else
             file_1 = fopen("base_thakkar.txt", "r");
           do {
             c = fscanf(file_1,"%s",oneword);
             if (strcmp(oneword,elemento) == 0 ){
               for (i = 1; i <=4 ; i++) 
                 d = fscanf(file_1,"%s",twoword); 
               fscanf(file_1,"%d",&entero1);
               printf("------ Basis set functions -------\n");
               printf("Internal basis set: %s\n", read_base);
               do {
                 d = fscanf(file_1,"%s",twoword);
                 if ( strcmp(twoword,"**") == 0 ) {
                   test_input = input_2_int(read_base, 
					    entero1, 
                                            Rc_temp,
                                            np, 
            				    mang, 
					    ncm, 
				   	    expo, 
					    &entero2, 
					    opt,
					    opt_flag,
					    file_1,
                                            nombre);
                   if (test_input == 1) return 1;
                   printf("Atomic Orbitals = %d\n", entero2);
                   printf("----------------------------------\n");
                   *base = entero2;
              //    if (strcmp(opt,"opt") == 0) j = 1;
              //    else j = 0;
              //    for (i = 0; i < entero2; i++) opt_flag[i] = j;
                   d = -1;
                 }
               } while (d == 1);
               c = -1;
             } 
           } while (c == 1 && c != EOF);
           fclose(file_1);
       } //End If checking Z
     } else {
         if (strcmp(read_base,"extern") == 0) {
           printf("Number of functions in the basis set: ");
           scanf("%d", &entero1);
           *base = entero1;
           printf("%d functions to read in the external basis set\n", *base);
           printf("------ Basis set functions -------\n");
           printf("External basis set\n");
           	  test_input = input_2_ext(read_base, 
					   entero1,
                                           Rc_temp,
                                           np, 
				   	   mang, 
					   ncm, 
					   expo, 
					   &entero2, 
					   opt_flag,
                                    	   file_1, 
				           nombre);
           if (test_input == 1) return 1;
           printf("Atomic Orbitals = %d\n", entero2);
           printf("----------------------------------\n");
           *base = entero2;
         } else {
             printf("You must write bunge, clementi, thakkar or extern to get a basis set functions\n");
             return 1;
         }
             
     } //End If, checking intern or extern
   }// End If, checking alpha and beta electrons
 return test_input;
 }


