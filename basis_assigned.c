#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int input_2_int(char   *read_base, 
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
		char   *nombre)
{
 int i, j, k, nt_temp;
 double temp;
 char config[4];

 FILE *archivo_tmp;

 
 time_t t;
 struct tm tstruct;
 char name_0[80];
 char name_1[80];
 char name_2[80];
 char name_3[80];
 int unsigned rand_1;
 srand((unsigned) time(NULL));

//mrb system("date | cut -c1-3,5-7,10-12,13-13,15-16,18-19,24-29 > date.dat");
//mrb strcpy(name_0,"date.dat");
//mrb archivo_tmp = fopen(name_0,"r");
//mrb fscanf(archivo_tmp, "%s %s %s", name_1, name_2, name_3);
//mrb fclose(archivo_tmp);
//mrb system("rm -f date.dat");

//mrb rand_1 = rand()%(3000 + 1);
//mrb sprintf(nombre,"tmp.%s_%s_%s_%d",name_1, name_2, name_3, rand_1);
 
 sprintf(nombre,"tmp.%lf", Rc);
 
//mrb  t = time(0);
//mrb  tstruct = *localtime(&t);
//mrb  strcpy(nombre,"tmp.");
//mrb  strcat(nombre,ctime(&t));
//
   archivo_tmp = fopen(nombre,"w");

// Reading basis set

       if (strcmp(opt,"opt") == 0) j = 1;
       else j = 0;


   nt_temp = -1;
   for (i = 0; i < base; i++) {
     fscanf(file_1,"%s %lf", config, &temp);
     printf("%s %12.6f\n", config, temp);
     fprintf(archivo_tmp,"%s %12.6f %d\n", config, temp, j);
     if (strcmp(config,"1S") == 0) {
       nt_temp++;
       np_temp[nt_temp] = 1;
       mang_temp[nt_temp] = 0;
       ncm_temp[nt_temp] = 0;
       expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
     } else if (strcmp(config,"2S") == 0) {
       nt_temp++;
       np_temp[nt_temp] = 2;
       mang_temp[nt_temp] = 0;
       ncm_temp[nt_temp] = 0;
       expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
     } else if (strcmp(config,"3S") == 0) {
       nt_temp++;
       np_temp[nt_temp] = 3;
       mang_temp[nt_temp] = 0;
       ncm_temp[nt_temp] = 0;
       expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
     } else if (strcmp(config,"4S") == 0) {
       nt_temp++;
       np_temp[nt_temp] = 4;
       mang_temp[nt_temp] = 0;
       ncm_temp[nt_temp] = 0;
       expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
     } else if (strcmp(config,"5S") == 0) {
       nt_temp++;
       np_temp[nt_temp] = 5;
       mang_temp[nt_temp] = 0;
       ncm_temp[nt_temp] = 0;
       expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
     } else if (strcmp(config,"6S") == 0) {
       nt_temp++;
       np_temp[nt_temp] = 6;
       mang_temp[nt_temp] = 0;
       ncm_temp[nt_temp] = 0;
       expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
     } else if (strcmp(config,"7S") == 0) {
       nt_temp++;
       np_temp[nt_temp] = 7;
       mang_temp[nt_temp] = 0;
       ncm_temp[nt_temp] = 0;
       expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
     } else if (strcmp(config,"2P") == 0) {
         for (k = 1; k <= 3; k++) {
           nt_temp++;
           np_temp[nt_temp] = 2;
           mang_temp[nt_temp] = 1;
           ncm_temp[nt_temp] = k - 2;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"3P") == 0) {
         for (k = 1; k <= 3; k++) {
           nt_temp++;
           np_temp[nt_temp] = 3;
           mang_temp[nt_temp] = 1;
           ncm_temp[nt_temp] = k - 2;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"4P") == 0) {
         for (k = 1; k <= 3; k++) {
           nt_temp++;
           np_temp[nt_temp] = 4;
           mang_temp[nt_temp] = 1;
           ncm_temp[nt_temp] = k - 2;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"5P") == 0) {
         for (k = 1; k <= 3; k++) {
           nt_temp++;
           np_temp[nt_temp] = 5;
           mang_temp[nt_temp] = 1;
           ncm_temp[nt_temp] = k - 2;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"6P") == 0) {
         for (k = 1; k <= 3; k++) {
           nt_temp++;
           np_temp[nt_temp] = 6;
           mang_temp[nt_temp] = 1;
           ncm_temp[nt_temp] = k - 2;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"7P") == 0) {
         for (k = 1; k <= 3; k++) {
           nt_temp++;
           np_temp[nt_temp] = 7;
           mang_temp[nt_temp] = 1;
           ncm_temp[nt_temp] = k - 2;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"3D") == 0) {
         for (k = 1; k <= 5; k++) {
           nt_temp++;
           np_temp[nt_temp] = 3;
           mang_temp[nt_temp] = 2;
           ncm_temp[nt_temp] = k - 3;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"4D") == 0) {
         for (k = 1; k <= 5; k++) {
           nt_temp++;
           np_temp[nt_temp] = 4;
           mang_temp[nt_temp] = 2;
           ncm_temp[nt_temp] = k - 3;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"5D") == 0) {
         for (k = 1; k <= 5; k++) {
           nt_temp++;
           np_temp[nt_temp] = 5;
           mang_temp[nt_temp] = 2;
           ncm_temp[nt_temp] = k - 3;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"6D") == 0) {
         for (k = 1; k <= 5; k++) {
           nt_temp++;
           np_temp[nt_temp] = 6;
           mang_temp[nt_temp] = 2;
           ncm_temp[nt_temp] = k - 3;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"7D") == 0) {
         for (k = 1; k <= 5; k++) {
           nt_temp++;
           np_temp[nt_temp] = 7;
           mang_temp[nt_temp] = 2;
           ncm_temp[nt_temp] = k - 3;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"4F") == 0) {
         for (k = 1; k <= 7; k++) {
           nt_temp++;
           np_temp[nt_temp] = 4;
           mang_temp[nt_temp] = 3;
           ncm_temp[nt_temp] = k - 4;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"5F") == 0) {
         for (k = 1; k <= 7; k++) {
           nt_temp++;
           np_temp[nt_temp] = 5;
           mang_temp[nt_temp] = 3;
           ncm_temp[nt_temp] = k - 4;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"6F") == 0) {
         for (k = 1; k <= 7; k++) {
           nt_temp++;
           np_temp[nt_temp] = 6;
           mang_temp[nt_temp] = 3;
           ncm_temp[nt_temp] = k - 4;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"7F") == 0) {
         for (k = 1; k <= 7; k++) {
           nt_temp++;
           np_temp[nt_temp] = 7;
           mang_temp[nt_temp] = 3;
           ncm_temp[nt_temp] = k - 4;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"5G") == 0) {
         for (k = 1; k <= 9; k++) {
           nt_temp++;
           np_temp[nt_temp] = 5;
           mang_temp[nt_temp] = 4;
           ncm_temp[nt_temp] = k - 5;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"6G") == 0) {
         for (k = 1; k <= 9; k++) {
           nt_temp++;
           np_temp[nt_temp] = 6;
           mang_temp[nt_temp] = 4;
           ncm_temp[nt_temp] = k - 5;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"7G") == 0) {
         for (k = 1; k <= 9; k++) {
           nt_temp++;
           np_temp[nt_temp] = 7;
           mang_temp[nt_temp] = 4;
           ncm_temp[nt_temp] = k - 5;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"6H") == 0) {
         for (k = 1; k <= 11; k++) {
           nt_temp++;
           np_temp[nt_temp] = 6;
           mang_temp[nt_temp] = 5;
           ncm_temp[nt_temp] = k - 6;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"7H") == 0) {
         for (k = 1; k <= 11; k++) {
           nt_temp++;
           np_temp[nt_temp] = 7;
           mang_temp[nt_temp] = 5;
           ncm_temp[nt_temp] = k - 6;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else if (strcmp(config,"7I") == 0) {
         for (k = 1; k <= 13; k++) {
           nt_temp++;
           np_temp[nt_temp] = 7;
           mang_temp[nt_temp] = 6;
           ncm_temp[nt_temp] = k - 7;
           expo_temp[nt_temp] = temp;
       opt_flag[nt_temp] = j;
         }
     } else {
          printf("Check basis set because I cannot use these orbitals\n");
          return 1;
       }
   }
   nt_temp++; // Total number of basis set functions

   *nt = nt_temp;

   fclose(archivo_tmp);
 return 0;
// End of basis set
}
