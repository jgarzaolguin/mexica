// Two-electron coulomb matrix
//
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//#ifdef _OPENMP
// #include <omp.h>
//#else
// #define omp_get_thread_num() 0
//#endif

void matrixg_coul(int nt, double* matp, double* matg, int* np,
                  int* mang, int* ncm, double* expo, double Rc, char* bound,
                  double* integrales_bie)
{
 int i, j, k, m, elemento1, elemento2, element_bie;

 int k2, cont, index_i, index_j, element_k, total_elements, proc_come;
 int bloque, restante, do_ini, do_fin, total_size, sizebi;
 double part_real;
 double *array_temp, *array_temp_2;

 double sum1, sum2, num1, num2;
 extern int indexes(int, int, int*, int*);
//mrb extern void memoria_double_uni(int, double**, char*);
 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);

 total_elements = nt*nt;

//#pragma omp parallel shared(nt, total_elements, integrales_bie, matg, matp) private(k2, index_i, index_j, i, j, k, m, num1, sum1, sum2, element_bie, elemento1)
//{ // begins omp
// int TID = omp_get_thread_num();
// #pragma omp for
 for (k2 = 0; k2 < total_elements; k2++) {
   indexes(nt, k2, &index_i, &index_j);
   index_i = index_i - 1;
   index_j = index_j - 1;
   i = index_i;
   j = index_j;
   sum1 = (double)0;
   for (k = 0; k < nt; k++) {
     sum2 = (double)0;
     for (m = 0; m < nt; m++)  {
       element_bie = i*nt*nt*nt + j*nt*nt + k*nt + m;
       num1 = integrales_bie[element_bie];
       elemento1 = m*nt + k;
       sum2 += matp[elemento1]*num1;
     }
     sum1 = sum1 + sum2;
   }
   matg[k2] = sum1;
  // printf("\n %f \n",matg[k2]);
 } 
//} // ends omp

 }

