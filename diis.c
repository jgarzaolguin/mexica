#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


int diis_check(int     nt, 
               double *matp,  
               double *mats,
               double *matfock, 
               double *matx, 
               double *mat_temp,
               double *diis_err)
{
 int i, valor, sizebi, dim;

 double zero, coef1, coef2, doble1, result, suma;
 double *mat_temp_1, *mat_temp_2, *mat_temp_3, *mat_resi;

 extern void summat(int, 
                    double, 
                    double*, 
                    double, 
                    double*, 
                    double*);
// extern void memoria_double_uni(int, 
//                                double**, 
//                                char*);
 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);

 dim = nt*nt;
 sizebi = dim*sizeof(double);

 valor = 0;
 coef1 = (double)1.0;
 coef2 = -(double)1.0;
 zero = (double)0.;

 mat_temp_1 =NULL;
 memoria_double_uni(sizebi, &mat_temp_1, "Mat_temp_1");
 mat_temp_2 =NULL;
 memoria_double_uni(sizebi, &mat_temp_2, "Mat_temp_2");
 mat_temp_3 =NULL;
 memoria_double_uni(sizebi, &mat_temp_3, "Mat_temp_3");
 mat_resi =NULL;
 memoria_double_uni(sizebi, &mat_resi, "Mat_resi");

 dgemm_("N", "N", &nt, &nt, &nt, &coef1, matp, &nt, mats, &nt, &zero,
        mat_temp_1, &nt);
 dgemm_("N", "N", &nt, &nt, &nt, &coef1, matfock, &nt, mat_temp_1, &nt, &zero,
        mat_temp_2, &nt);

 dgemm_("N", "N", &nt, &nt, &nt, &coef1, matp, &nt, matfock, &nt, &zero,
        mat_temp_1, &nt);
 dgemm_("N", "N", &nt, &nt, &nt, &coef1, mats, &nt, mat_temp_1, &nt, &zero,
        mat_temp_3, &nt);

 summat(dim, coef1, mat_temp_2, coef2, mat_temp_3, mat_resi);

 dgemm_("N", "N", &nt, &nt, &nt, &coef1, mat_resi, &nt, matx, &nt, &zero,
        mat_temp_1, &nt);
 dgemm_("T", "N", &nt, &nt, &nt, &coef1, matx, &nt, mat_temp_1, &nt, &zero,
        mat_temp_2, &nt);

 result = -(double)1.0;
 for ( i = 0; i < nt*nt; i++) {
   doble1 = fabs(mat_temp_2[i]);
   result = ( doble1 > result ? doble1 : result );
 }

 for ( i = 0; i < nt*nt; i++) mat_temp[i] = mat_temp_2[i];

 *diis_err = result;


 free(mat_resi);
 free(mat_temp_3);
 free(mat_temp_2);
 free(mat_temp_1);

 return (valor);
 }

int solving_C(int new_dim, double *mat_B, double *vect_1, double *X)
{
 int i, j, k, new_size, element2;
 int NRHS, LDB, LDX, INFO;
 int *IPIV,*IWORK;
 double RCOND;
 double *AP, *AFP, *FERR, *BERR, *WORK;

// extern void memoria_double_uni(int, 
//                                double**, 
//                                char*);
 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);

 NRHS = 1;
 new_size = ((new_dim*(new_dim + 1))/2)*sizeof(double);
 AP =NULL;
 memoria_double_uni(new_size, &AP, "AP");

 k = -1;
 for ( j = 0; j < new_dim; j++)
   for ( i = 0; i < new_dim; i++)
     if (j >= i) {
       k++;
       element2 = i*new_dim + j;
       AP[k] = mat_B[element2];
     }

 AFP =NULL;
 memoria_double_uni(new_size, &AFP, "AP");
 new_size = new_dim*sizeof(int);
 IPIV =NULL;
 memoria_integer_uni(new_size, &IPIV, "IPIV");
 IWORK =NULL;
 memoria_integer_uni(new_size, &IWORK, "IWORK");
 LDB = new_dim;
 new_size = new_dim*sizeof(double);
 LDX = new_dim;
 new_size = NRHS*sizeof(double);
 FERR =NULL;
 memoria_double_uni(new_size, &FERR, "FERR");
 BERR =NULL;
 memoria_double_uni(new_size, &BERR, "BERR");
 new_size = 3*new_dim*sizeof(double);
 WORK =NULL;
 memoria_double_uni(new_size, &WORK, "WORK");

 dspsvx_ ("N","U", &new_dim, &NRHS, AP, AFP, IPIV, vect_1, &LDB, X,
        &LDX, &RCOND, FERR, BERR, WORK, IWORK, &INFO );

 free(AP);
 free(AFP);
 free(IPIV);
 free(IWORK);
 free(FERR);
 free(BERR);
 free(WORK);

 return (INFO);

 }

//
//

int main_diis(int     nt, 
              int     iter_inter, 
              double *mat_temp_1, 
              double *mat_temp_2,
              double *mat_temp_3, 
              double *mat_e_store[100], 
              double *mat_f_store[100],
              double *matfock)
{
 int i, j, m, k, element2, new_dim, new_size;
 double suma, coef1, zero;
 double *mat_B, *vect_1, *vect_c;

// extern void memoria_double_uni(int, 
//                                double**, 
//                                char*);
 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);
 extern int solving_C(int, 
                      double*, 
                      double*, 
                      double*);

 coef1 = (double)1.;
 zero = (double)0.;
 new_dim = iter_inter + 1;
 new_size = (new_dim*new_dim)*sizeof(double);
 mat_B =NULL;
 memoria_double_uni(new_size, &mat_B, "Mat_B");
 new_size = new_dim*sizeof(double);
 vect_1 =NULL;
 memoria_double_uni(new_size, &vect_1, "vect_1");
// Building mat_B and vector_C
 vect_1[0] = -(double)1.;
 mat_B[0] = (double)0.;
 for (i = 1; i <= iter_inter; i++)
   for (j = 1; j <= iter_inter; j++) {
     for (m = 0; m < nt*nt; m++) {
       mat_temp_1[m] = mat_e_store[i][m];
       mat_temp_2[m] = mat_e_store[j][m];
     }
     dgemm_ ("N", "T", &nt, &nt, &nt, &coef1, mat_temp_1, &nt, mat_temp_2, &nt, &zero,
            mat_temp_3, &nt);
     suma = (double)0.;
     for (m = 0; m < nt; m++) suma = suma + mat_temp_3[m*nt + m];
     element2 = i*new_dim + j;
     mat_B[element2] = suma;
   }

  for (i = 1; i <= iter_inter; i++) {
     vect_1[i] = (double)0.;
     mat_B[i] = -(double)1.0;
     mat_B[i*(iter_inter + 1)] = -(double)1.0;
  }

 new_size = new_dim*sizeof(double);
             vect_c = NULL;
             memoria_double_uni(new_size, &vect_c, "vect_c");

             solving_C(new_dim, mat_B, vect_1, vect_c);

             for ( i = 0; i < nt*nt; i++) matfock[i] = (double)0.;

             for (k = 1; k <= iter_inter; k++)
               for ( i = 0; i < nt*nt; i++)
                 matfock[i] = matfock[i] + vect_c[k]*mat_f_store[k][i];

             free(vect_c);
             vect_c = 0;
             free(vect_1);
             vect_1 = 0;
             free(mat_B);
 return (0);
}

