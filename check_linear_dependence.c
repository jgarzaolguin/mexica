#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

 int checking_linear_dependece(int       nt,
                               int      *bt,
                               double   *expo,
                               int      *mang,
                               char     *nombre)
 {//label 1
 
  int deg,
      ang_test,
      bandera, 
      index,
      todos,
      i,
      j,
      k;
  int    result;
  double scratch,
         diff;
  int    scratch2;
  char   config[4];
  char   symbol[200][4];
  int    save_deg[200];
  FILE*  archivo_opt;

     i =-1;
     j = 0;
     result = 0;

     archivo_opt = fopen(nombre,"r");

     while (fscanf(archivo_opt,"%s %lf %d", config, &scratch, &scratch2) != EOF)
     {
      strcpy(symbol[j], config);
      if (config[1] == 'S') i = i + 1;
      if (config[1] == 'P') i = i + 3;
      if (config[1] == 'D') i = i + 5;
      if (config[1] == 'F') i = i + 7;
      if (config[1] == 'G') i = i + 9;
      if (config[1] == 'H') i = i + 11;
      if (config[1] == 'I') i = i + 13;
       save_deg[j] = i;
       //printf("MRB  %s  >>%s  %lf >>%d %d\n", &config[0], symbol[j], expo[i], save_deg[j], j);
     j++;
     }
     fclose(archivo_opt);
    
     *bt = j;

     for (i = 0; i < j - 1; i++) 
      for (k = i + 1; k < j; k++) {//label 2 
         
         if (expo[save_deg[i]] < 0.01 || 
             expo[save_deg[k]] < 0.01 || 
             expo[save_deg[i]] < 0.f || 
             expo[save_deg[k]] < 0.f ){
            printf("\nvalue exponent down to 0.1 >>  %s %lf %d>> %s %lf %d\n", symbol[i], 
                                                                               expo[save_deg[i]], 
                                                                               i, 
                                                                               symbol[k], 
                                                                               expo[save_deg[k]], 
                                                                               k);
            result = 1;
            break;
           } 
         if (strcmp(symbol[i], symbol[k]) == 0) {
           diff = fabs(expo[save_deg[i]]-expo[save_deg[k]]);
           if (diff < 0.01) {
              printf("\ncheck_dependece difference >>  %lf\n", diff);
              printf("check_dependence  >>%s  %lf >>%d %d\n", symbol[i], expo[save_deg[i]], save_deg[i], i);
              printf("check_dependence  >>%s  %lf >>%d %d\n\n", symbol[k], expo[save_deg[k]], save_deg[k], k);
              expo[save_deg[i]] = expo[save_deg[i]];
              expo[save_deg[k]] = expo[save_deg[k]];
              result = 1;
              break;
             }

           }
        }//label 2

 return (result);

 }//label 1
                              
                              
                              
