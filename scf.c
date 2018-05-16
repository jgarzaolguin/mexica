// Distributing two-electron integrals over OpenMP
// Computing total energy
// SCF process
//
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#ifdef _OPENMP
 #include <omp.h>
#else
 #define omp_get_thread_num() 0
#endif

//
// Function to distribute two-electron integrals over several cores
// by using OpenMP programming techniques
//
extern  int bielectronicas_CPU(char   *using_gamma,
                               int     nt, 
                               int    *np, 
                               int    *mang, 
                               int    *ncm, 
                               double *expo,
                               char   *bound, 
                               double  Rc, 
                               double  gamma_couple, 
                               double *integrales_bie, 
                               double *zeta,
                               double *N_minus, 
                               double *N_plus,
                               double *arreglo_factorial, 
                               double *arreglo_inv_factorial)
{
  int i, indice_i, indice_j, indice_k, indice_l;
  int cuad, cubo, cuart, high_l;
  double integral_cpu, inv_Rc, pi, *two_l_plus_1;

  extern double doselec(char *using_gamma, int mu, int nu, int lam, int sig, double Rc,
                        double inv_Rc, double gamma_couple, char *bound, double *expo, int *np,
                        int *ang, int    *ncm, double *zeta, double *N_minus, double *N_plus,
                        double *arreglo_factorial, double *arreglo_inv_factorial,
                        double *two_l_plus_1);

  pi = 4.f*atan(1.f);
  high_l = 50;
  two_l_plus_1 = (double *)malloc(high_l*sizeof(double));

  for (i = 0; i < high_l; i++)
    two_l_plus_1[i] = 4.f*pi/((double) 2*i + 1.f);
  cuad = nt*nt;
  cubo = cuad*nt;
  cuart = cubo*nt;
  inv_Rc = 1.f/Rc;
  #pragma omp parallel shared(integrales_bie, Rc, bound, expo, np, mang, ncm, cuad, cubo, cuart, nt, gamma_couple, two_l_plus_1) private(i, indice_i,indice_j, indice_k, indice_l, integral_cpu)
  {
    #pragma omp for
    for (i = 0; i < cuart; i++) {
       indice_i = i/cubo;
       indice_j = i - indice_i*cubo;
       indice_j = indice_j/cuad;
       indice_k = i - indice_i*cubo - indice_j*cuad;
       indice_k = indice_k/nt;
       indice_l = i - indice_i*cubo - indice_j*cuad - indice_k*nt;       

       integral_cpu = doselec(using_gamma,
                              indice_i, 
                              indice_j, 
                              indice_k, 
                              indice_l, 
                              Rc, 
                              inv_Rc,                               
                              gamma_couple,  
                              bound, 
                              expo, 
                              np, 
                              mang, 
                              ncm, 
                              zeta, 
                              N_minus, 
                              N_plus, 
                              arreglo_factorial, 
                              arreglo_inv_factorial, two_l_plus_1);

       integrales_bie[i] = integral_cpu;                
     }                                                                                                 
    }  // Termina openmp 
  free(two_l_plus_1);
  return 0;
}
//
// Function to compute total energy
//
int compute_energy(int     compara, 
		   int     nt, 
		   double *matp, 
		   double *matk,
		   double *matkint,
		   double *matkext,
                   double *matv, 
		   double *matg_coul, 
		   double *matg_exch,
                   double *matpalfa, 
		   double *matpbeta,
                   double *matg_exch_alfa, 
		   double *matg_exch_beta,
                   double *e_kin, 
                   double *e_kinint, 
                   double *e_kinext, 
		   double *e_v, 
		   double *e_core, 
		   double *e_coul,
                   double *e_exch, 
		   double *energia)
{
 int i, j, element1, element2;
 double sumak, sumakint, sumakext, sumav, suma2, suma3, doble1;

 sumak    = (double)0.;
 sumakint = (double)0.;
 sumakext = (double)0.;
 sumav    = (double)0.;
 suma2    = (double)0.;
 suma3    = (double)0.;

 for (i = 0; i < nt; i++)
   for (j = 0; j < nt; j++) {
     element1 = j*nt + i;
     element2 = i*nt + j;
     sumak    += matp[element2]*matk[element1];
     sumakint += matp[element2]*matkint[element1];
     sumakext += matp[element2]*matkext[element1];
     sumav    += matp[element2]*matv[element1];
     suma2    += (double)0.5*matp[element2]*matg_coul[element1];
   }

 if (compara == 0) {
   for (i = 0; i < nt; i++)
     for (j = 0; j < nt; j++) {
       element1 = j*nt + i;
       element2 = i*nt + j;
       suma3 += (double)0.5*matp[element2]*matg_exch[element1];
     }
 } else {
   for (i = 0; i < nt; i++)
     for (j = 0; j < nt; j++) {
       element1 = j*nt + i;
       element2 = i*nt + j;
       doble1 = matpalfa[element2]*matg_exch_alfa[element1];
       suma3 += (double) 0.5*(doble1 + matpbeta[element2]*matg_exch_beta[element1]);
     }
 }
 *e_kin    = sumak;
 *e_kinint = sumakint;
 *e_kinext = sumakext;
 *e_v      = sumav;
 *e_core   = sumak + sumav;
 *e_coul   = suma2;
 *e_exch   = suma3;
 *energia  = (sumak + sumav + suma2);

  return(0);
 }
//
///////////
//
// Function to control the SCF process
//
extern  int scf(int     nt, 
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
                double *total_energy, 
                int     print_vectors, 
                double  epsilon,
                int     imprime,
                int     plasma,
                double *cusp_kato,
                char *properties)
{
 float          tiempo_usado;
 int i, j, selected_orb, compara, iter, iter_inter, dim, elec, valor, bandera;
 int sizedouble, sizebi; //bandera_GPU ;
 int abre;
 int h;
 unsigned long long int cuad, cubo, cuart;
 double arreglo_factorial[100];
 double arreglo_inv_factorial[100];
double rho_int,  
       rho,
       rho_ini,
       der_rho_int,
       div_rho_int,
       rdf_int;
double rho_ext, 
       der_rho_ext,
       div_rho_ext,
       rdf_ext;
double shannon,
       shannon_i,
       radius_i,
       slope,
       intercept,
       bound_rho,
       shannon_entropy;

 int n_points;
 double step_r, radius, radiuscons, numero_pi;

 double doble1;

      doble1 = (double)1;
      numero_pi = atan(doble1)*(double)4;


 FILE *write_out;
 char name[50];

 cuad = (unsigned long long int) nt*nt;
 cubo = (unsigned long long int) nt*cuad;
 cuart = (unsigned long long int) nt*cubo;
 double coef1, 
        e_kin, e_kinint, e_kinext, e_v, e_check;
 double energia, energiavieja, e_core, e_coul, e_exch, difer, difmezcla;
 double max_e, max_e_beta, max_e_tope, mp2_energy, omega;
 double *integrales_bie;
 double *arreglo_mo,
        *arreglo_mo_ab;
 //mrb double *integrales_post_HF;
 //mrb int    *integral_importante;
 int    integrales;
 double tolbie;

 double *valores, *valoresalfa, *valoresbeta;
 double *vectsfin, *vectsfinalfa, *vectsfinbeta, *vect_alfa_post, *vect_beta_post;
 double *mats, *matx, *matk, *matkint, *matkext, *matv, *matcore;
 double *matp, *matpalfa, *matpbeta;
 double *matgvieja, *matg, *matg_alfa, *matg_beta, *matgvieja_alfa, *matgvieja_beta;
 double *matg_coul, *matg_exch, *matg_exch_alfa, *matg_exch_beta;
 double *matfock, *matfockalfa, *matfockbeta, *mat_temp_1, *mat_temp_2, *mat_temp_3;
 double *mat_temp_1_beta, *mat_temp_2_beta, *mat_temp_3_beta;
 double *mat_f_store[800], *mat_e_store[800];
 double *mat_f_store_beta[800], *mat_e_store_beta[800];
 double *NC_minus, *NC_plus, *zetas, tmp1, tmp2, tmp3;
 double energia_prueba;
 double SHAN;
 int    n_points_temp;


 extern	void traslape(char   *using_gamma,
                      int     nt,
	              double *mats, 
	              int    *np, 
	              int    *mang,
	              int    *ncm, 
	              double *expo,
                      char   *basis, 
	              double  Rc,
	              double  gamma_couple,
	              char   *bound, 
	              double  U_0,
	              double *NC_minus, 
	              double *NC_plus, 
	              double *arreglo_factorial,
	              double *arreglo_inv_factorial);


//mrb extern void transforma(int, 
//mrb                        double* , 
//mrb                        double* 
//mrb                        int*);
 extern void transforma(int     nt,
                        double *mats, 
                        double *matx);


//mrb extern void valorespropios(int, 
//mrb                            double*, 
//mrb                            double*, 
//mrb                            double*);
 extern void valorespropios(int     nt,
                            double *mata,
                            double *valores,
                            double *vectores,
                            int    *revisa);



 extern	void cinetica(char   *using_gamma,
                      int     nt,
	              double *matk, 
	              double *matkint, 
	              double *matkext, 
                      char   *basis,
	              int    *np, 
	              int    *mang, 
	              int    *ncm, 
	              double *expo, 
	              double  Rc,
	              double  gamma_couple,
	              char   *bound,
	              double *NC_minus, 
	              double *NC_plus, 
	              double *arreglo_factorial,
	              double *arreglo_inv_factorial);


 extern	void potencial(char   *using_gamma,
                       int     nt,
	               double  z,
	               double *matv, 
	               int    *np, 
	               int    *mang,
	               int    *ncm, 
	               double *expo, 
                       char   *basis,
	               double  Rc,
	               double  gamma_couple,
	               char   *bound, 
                       int     iter_pol,
                       double  charge_int,
	               double  U_0,
	               double *NC_minus, 
	               double *NC_plus, 
	               double *arreglo_factorial,
	               double *arreglo_inv_factorial,
                       double elec,
                       int plasma);


 extern void newfock(int, 
                     double*, 
                     double*, 
                     double*);

 extern void matdens(int, 
                     int, 
                     double*, 
                     double*, 
                     char*);

 extern void matrixg_coul(int, 
                          double*, 
                          double*, 
                          int*, 
                          int*, 
                          int*,  
                          double*,
                          double, 
                          char*, 
                          double*);

 extern void matrixg_exch(int, 
                          double*, 
                          double*, 
                          int*, 
                          int*, 
                          int*,
                          double*, 
                          double, 
                          char*, 
                          char*, 
                          double*);

//mrb extern void eigensystem_matdens(int , 
//mrb                                 char*, 
//mrb                                 int, 
//mrb                                 double*, 
//mrb                                 double*,
//mrb                                 double*, 
//mrb                                 double*, 
//mrb                                 double*);

 extern void eigensystem_matdens(int     nt,
                                char   *tipo,
                                int     elec,
                                double *matx,
                                double *matfock,
                                double *vals,
                                double *vectsfin,
                                double *matp,
                                int    *revisa_1);



 extern int diis_check(int     nt,
                       double *matp,
                       double *mats,
                       double *matfock,
                       double *matx,
                       double *mat_temp,
                       double *diis_err);


 extern int solving_C(int, 
                      double*, 
                      double*, 
                      double*);

/*
 extern int compute_energy(int, int, double*, double*, double*, double*, double*,
                           double*, double*, double*, double*, double*, double*,
                           double*, double*, double*, double*);
*/
extern int compute_energy(int     compara,
                          int     nt,
                          double *matp,
                          double *matk,
                          double *matkint,
                          double *matkext,
                          double *matv,
                          double *matg_coul,
                          double *matg_exch,
                          double *matpalfa,
                          double *matpbeta,
                          double *matg_exch_alfa,
                          double *matg_exch_beta,
                          double *e_kin,
                          double *e_kinint,
                          double *e_kinext,
                          double *e_v,
                          double *e_core,
                          double *e_coul,
                          double *e_exch,
                          double *energia);


 extern int exact_exch_ener(int, 
                            int, 
                            double*, 
                            double*, 
                            double* , 
                            double* ,
                            double* , 
                            double* , 
                            double* );

 extern int main_diis(int     nt,
                      int     iter_inter,
                      double *mat_temp_1,
                      double *mat_temp_2,
                      double *mat_temp_3,
                      double *mat_e_store[100],
                      double *mat_f_store[100],
                      double *matfock);

// extern int rho_derrho_on_mesh(int nt, double Rc, double* matp, int bandera_GPU,
//                        double* expo, int* np, int* ang, int* ncm,
//                        double infinity, int espacios,
//                        double* rho, double* derrho);

 extern int bielectronicas_GPU(int, 
                               int*, 
                               int*, 
                               int*,
                               double*, 
                               double,
                               char*, 
                               int, 
                               double*);

extern int mp2_cpu(int nt, 
                   int compara, 
                   int elecalfa, 
                   int elecbeta,
                   double *vectalfa, 
                   double *vectbeta,
                   double *valalfa, 
                   double *valbeta, 
                   double *integrales_bie,
                   double *arreglo_mo,
                   double *total_energy);

extern int ep2_cpu(char   *espin,
                   int     nt,
                   int     elecalfa,
                   int     elecbeta,
                   double *arreglo,
                   double *vectalfa,
                   double *vectbeta,
                   double *valalfa,
                   double *valbeta,
                   int     orbital,
                   double *arreglo_mo,
                   double *arreglo_mo_ab,
                   double *total_energy);


 extern void expected_value_r(int     nt,
                              int     r_exp,
                              int    *np,
                              int    *mang,
                              int    *ncm,
                              double *expo,
                              double  Rc,
                              char   *bound,
                              char   *basis,
                              double *vectsfin,
                              int     elecalfa,
                              double  gamma_couple,
                              double *NC_minus,
                              double *NC_plus,
                              double *arreglo_factorial,
                              double *arreglo_inv_factorial,
                              char *using_gamma,
                              double *zeta,
                              double *grid);




 extern double confined_orbital(int     nt, 
                                int     orbital, 
                                double  r, 
                                double  Rc, 
                                double *expo,
                                int    *np, 
                                double *vectors);

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

 extern int shannon_entropy_finite(int     z,
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
                                   double  gamma_couple);


  extern int shannon_entropy_confined(int     z,
                                      int     compara,
                                      int     nt,
                                      int     elecalfa,
                                      int     elecbeta,
                                      double  Rc,
                                      double *expo,
                                      int    *np,
                                      double *vectsfinalfa,
                                      double *vectsfinbeta,
                                      char   *tipo,
                                      double *arreglo_factorial,
                                      double *arreglo_inv_factorial);

 extern  void  print_spherical_density_finite(int     compara,
                                              char   *using_gamma,
                                              int     nt,
                                              double  z,
                                              int     elecalfa,
                                              int     elecbeta,
                                              double  Rc,
                                              double *expo,
                                              double *zetas,
                                              int    *np,
                                              int    *mang,
                                              double *vectsfinalfa,
                                              double *vectsfinbeta,
                                              char   *tipo,
                                              double  gamma_couple,
                                              double *NC_minus,
                                              double *NC_plus);


 extern  void  print_spherical_density_confined(int     compara,
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
                                                double *arreglo_inv_factorial);






 extern double rho_radial_free(int     nt, 
                               int     elec, 
                               double  r, 
                               double *expo,
                               int    *np, 
                               double *vectors, 
                               char   *tipo);

 extern double der_rho_radial_free(int     nt, 
                                   int     elec, 
                                   double  r, 
                                   double *expo,
                                   int    *np, 
                                   double *vectors, 
                                   char   *tipo);

 extern int constants_normalization_finite(int     mu, 
					   int    *np, 
			                   int    *ang, 
					   double *zetas, 
		                           double  Rc,
                                           double *arreglo_factorial, 
				           double *arreglo_inv_factorial,
                                           char   *using_gamma,
                             		   double  gamma_couple,
                             		   double *const_n_minus, 
			                   double *const_n_plus, 
			                   double *alphas, int print_expo);

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



 extern double  der_rho_radial_finite_int(char   *using_gamma,
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

 extern void summat(int, 
                    double, 
                    double*, 
                    double, 
                    double*, 
                    double*);

 extern void copiamat(int, 
                      double*, 
                      double*);

 extern void multmat(int, 
                     double*, 
                     double*, 
                     double*);

 extern void imprimemat(int     nt, 
                        double *mat);

//mrb extern void memoria_double_uni(int, 
//mrb                                double**, 
//mrb                                char*);

 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);

 extern int bielectronicas_CPU(char   *using_gamma,
                               int     nt, 
                               int    *np, 
                               int    *mang, 
                               int    *ncm, 
                               double *expo, 
                               char   *bound,
                               double  Rc, 
                               double  gamma_couple, 
                               double *integrales_bie,
                               double *zeta, 
                               double *N_minus, 
                               double *N_plus,
                               double *arreglo_factorial, 
                               double *arreglo_inv_factorial);
 extern double factorial(int);

 extern int grid_rhorad(int           z,
                        char         *using_gamma,
                        int           compara,
                        char         *bound,
                        char         *basis, /* new */
                        int           nt,
                        int           elecalfa,
                        int           elecbeta,
                        double        Rc,
                        double       *expo,
                        int          *np,
                        double       *zetas,
                        int          *mang,
                        double       *vectsfinalfa,
                        double       *vectsfinbeta,
                        char         *tipo,
                        double       *NC_minus,
                        double       *NC_plus,
                        double        gamma_couple,
                        double       *grid,
                        double       *grid_rho,
                        double       *grid_rho_beta,
                        int           n_points,
                        int           n_boundary,
                        double       *arreglo_factorial,
                        double       *arreglo_inv_factorial);

 extern int grid_rhorad_orbital(int           z,
                        char         *using_gamma,
                        int           compara,
                        char         *bound,
                        int           nt,
                        int           elecalfa,
                        int           elecbeta,
                        double        Rc,
                        double       *expo,
                        int          *np,
                        double       *zetas,
                        int          *mang,
                        double       *vectsfinalfa,
                        double       *vectsfinbeta,
                        char         *tipo,
                        double       *NC_minus,
                        double       *NC_plus,
                        double        gamma_couple,
                        double       *grid,
                        double       *grid_rho,
                        double       *grid_rho_beta,
                        int           n_points,
                        int           n_boundary,
                        double       *arreglo_factorial,
                        double       *arreglo_inv_factorial,
                        int           selected_orb);

 extern int grid_derrad(int           z,                                    
                        char         *using_gamma,
                        int           compara,
                        char         *bound,
                        char         *basis,
                        int           nt,
                        int           elecalfa,
                        int           elecbeta,
                        double        Rc,
                        double       *expo,
                        int          *np,
                        double       *zetas,
                        int          *mang,
                        double       *vectsfinalfa,
                        double       *vectsfinbeta,
                        char         *tipo,
                        double       *NC_minus,
                        double       *NC_plus,
                        double        gamma_couple,
                        double       *grid,
                        double       *grid_der,
                        double       *grid_der_beta,
                        int           n_points,
                        int           n_boundary,
                        double       *arreglo_factorial,
                        double       *arreglo_inv_factorial);

 extern int grid_secderrad(int           z,                                    
                           char         *using_gamma,
                           int           compara,
                           char         *bound,
                           int           nt,
                           int           elecalfa,
                           int           elecbeta,
                           double        Rc,
                           double       *expo,
                           int          *np,
                           double       *zetas,
                           int          *mang,
                           double       *vectsfinalfa,
                           double       *vectsfinbeta,
                           char         *tipo,
                           double       *NC_minus,
                           double       *NC_plus,
                           double        gamma_couple,
                           double       *grid,
                           double       *grid_secder,
                           double       *grid_secder_beta,
                           int           n_points,
                           int           n_boundary,
                           double       *arreglo_factorial,
                           double       *arreglo_inv_factorial);

 extern double xc_energy(double *correlationc,
                         double *exchangex,
                         int     points,
                         int     compara,
                         char   **save_dft,
                         double *weight_dft,
                         int     flag_dft,
                         double *energy_array,
                         double *grid,
                         double *rho_alpha,
                         double *rho_beta,
                         double *derho_alpha,
                         double *derho_beta,
                         double *d2rho_alpha,
                         double *d2rho_beta);
  
 extern void xc_potential (int           points,
                           int           compara,
                           int           point_interm,
                           char        **save_dft,
                           double       *weight_dft,
                           int           flag_dft,
                           char         *using_gamma,
                           char         *bound,
                           double        Rc,
                           int           nt,
                           double       *expo,
                           int          *np,
                           double       *zetas,
                           int          *mang,
                           int          *ncm,
                           double       *NC_minus,
                           double       *NC_plus,
                           double        gamma_couple,
                           double       *pot_array,
                           double       *pot_array_beta,
                           double       *local_int,
                           double       *local_int_beta,
                           double       *grid,
                           double       *rho_alpha,
                           double       *rho_beta,
                           double       *derho_alpha,
                           double       *derho_beta,
                           double       *d2rho_alpha,
                           double       *d2rho_beta,
                           double       *mat_ks,
                           double       *mat_ks_beta,
                           double       *arreglo_factorial,
                           double       *arreglo_inv_factorial,
                           double       *nolocal_pot_array);

 extern int wf_closed_shell(int z, char *using_gamma, int compara, char *bound, int nt, int elecalfa, int elecbeta,
                            double Rc, double *expo, int *np, double *zetas, int *mang, int *ncm, double *vectsfinalfa,
                            double *vectsfinbeta, char *tipo, double *NC_minus, double *NC_plus, double gamma_couple,
                            double *grid, double *grid_rho, double *grid_rho_beta, double *grid_der, double *grid_der_beta,
                            int n_points, 
                            double *arreglo_factorial, double *arreglo_inv_factorial,
                            double *matp, double *mats, int *iter, int save_i, int print_vectors, double *cusp_kato);



 extern double numerical_int(double   *grid,
                             double   *grid_fun,
                             int       points);
 extern int building_grid(int      z,
                          char    *bound,
                          double  Rc,
                          double *grid,
                          int     n_points,
                          int    *save_i);

 extern void xc_over_grid(int compara, char **save_dft, int flag_dft, double *weight_dft, int n_points,
                          double * grid, double *grid_rho, double *grid_der);

 extern int Evaluate_Elect_Pot(double z, int nt, double *matp, int *np, int *ang, int *ncm,
                      double *expo, char *bound, double *arreglo_factorial,
                      double *arreglo_inv_factorial, double *grid, int n_points, double Rc,
                      double *NC_minus, double *NC_plus, char *basis);

 extern int print_out_array(int points, double *grid, double *array, char *name_file);


 time_t time_scf_ini, time_scf_fin, time_bie_ini, time_bie_fin, time_3, time_4;
 double cx_D, 
        e_xc,
        temp_xcmat,
        coef_HF, 
        coef13, 
        coef43, 
        suma_x, 
        point_r, 
        *mat_ks,
        *mat_ks_beta,
        *grid, 
        *grid_rho, 
        *grid_rho_beta, 
        *grid_der, 
        *grid_secder, 
        *grid_secder_beta, 
        *grid_der_beta, 
        *array_i, *array_ii, *array_iii, *array_iv, *nolocal_pot_array;

 double Pi = 4.f*atan(1.f);

 int save_i;
 int save_coef_HF;
 coef13 = 1.f/3.f;
 coef43 = 4.f*coef13;
 cx_D = -0.7386f;
 coef_HF = 0.f;
 save_coef_HF = 0;
 n_points = 600;               //número de puntos
 grid = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &grid, "Grid");
 grid[0] = 0.f;

 double correlationc[flag_dft];
 double exchangex[flag_dft];

 for (i = 0; i < flag_dft; i++) {
 correlationc[i] = 0.f;
 exchangex[i]    = 0.f;
 }
//identifica Hartree Fock
 for (i = 1; i < flag_dft; i++) {
    if(strcmp(save_dft[i], "hf") == 0) {
      save_coef_HF = i;
    }
 }
 
 if (save_coef_HF != 0)
   coef_HF = weight_dft[save_coef_HF];
 else
   coef_HF = 0.f;

 save_i = 0;

 building_grid(z,
              bound,
              Rc,
              grid,
              n_points,
              &save_i);

 if(strcmp(bound,"confined") == 0)  n_points = save_i + 1;

 grid_rho = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &grid_rho, "Grid_rho");

 grid_rho_beta = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &grid_rho_beta, "Grid_rho_beta");

 grid_der = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &grid_der, "Grid_der");

 grid_der_beta = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &grid_der_beta, "Grid_der_beta");

 grid_secder = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &grid_secder, "Grid_sec_der");

 grid_secder_beta = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &grid_secder_beta, "Grid_sec_der_beta");

 array_i = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &array_i, "Array_i");

 array_ii = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &array_ii, "Array_ii");

 array_iii = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &array_iii, "Array_iii");

 array_iv = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &array_iv, "Array_iv");
 
 nolocal_pot_array = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &nolocal_pot_array, "nolocal_pot_array");

 char name_file_vectors[25];
 char name_file_vectors_beta[25];
 FILE* vectors_file;
 FILE* vectors_file_beta;

 time_3 = time (NULL);
 time_scf_ini = time (NULL);
 valor = 0;     // mike

 if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"rks") == 0) 
   compara = 0;
 else
   compara = 1;
 if (compara == 0) printf("Restricted wave-function\n");
 else              printf("Unrestricted wave function\n");

 if (compara == 0) {
 printf("\n DFT ");
 for (i = 1; i < flag_dft; i++)
    printf("  %s    %5.2lf", save_dft[i], weight_dft[i]);
 printf("\n\n");
 } else
    printf("  %s    %5.2lf", save_dft[i], weight_dft[i]);

 elec = elecalfa + elecbeta;

 dim = nt*nt;
 sizedouble = nt*sizeof(double);
 sizebi = dim*sizeof(double);

// /////////////////////////////////////////////////////////////////////
// Memory allocation for the arrays

// Only for processor 0
// if (myproc == 0) {
   mats =NULL;
   memoria_double_uni(sizebi, &mats, "Mats");
// One-electron matrix
   matk =NULL;
   memoria_double_uni(sizebi, &matk, "Cinetica");
   matkint =NULL;
   memoria_double_uni(sizebi, &matkint, "Cinetica interna");
   matkext =NULL;
   memoria_double_uni(sizebi, &matkext, "Cinetica externa");
   matv =NULL;
   memoria_double_uni(sizebi, &matv, "Matv");
   matcore =NULL;
   memoria_double_uni(sizebi, &matcore, "Matcore");
   for (i = 0; i < nt*nt; i++) {
      mats[i] = 0.f;
      matk[i] = 0.f;
      matv[i] = 0.f;
      matcore[i] = 0.f;
   }
// Two-electron matrix
   matg_coul =NULL;
   memoria_double_uni(sizebi, &matg_coul, "Matg_coul");
   for (i = 0; i < nt*nt; i++) matg_coul[i] = 0.f;
   if (compara == 0) {
     matgvieja =NULL;
     memoria_double_uni(sizebi, &matgvieja, "Matgvieja");
     matg =NULL;
     memoria_double_uni(sizebi, &matg, "Matg");
     matg_exch =NULL;
     memoria_double_uni(sizebi, &matg_exch, "Matg_exch");
     for (i = 0; i < nt*nt; i++) { 
      matgvieja[i] = 0.f;
      matg[i] = 0.f;
      matg_exch[i] = 0.f;
     }
   }
   else {
     matgvieja_alfa =NULL;
     memoria_double_uni(sizebi, &matgvieja_alfa, "Matgvieja_alfa");
     matgvieja_beta =NULL;
     memoria_double_uni(sizebi, &matgvieja_beta, "Matgvieja_beta");
     matg_alfa =NULL;
     memoria_double_uni(sizebi, &matg_alfa, "Matg_alfa");
     matg_beta =NULL;
     memoria_double_uni(sizebi, &matg_beta, "Matg_beta");
     matg_exch_alfa =NULL;
     memoria_double_uni(sizebi, &matg_exch_alfa, "Matg_exch");
     matg_exch_beta =NULL;
     memoria_double_uni(sizebi, &matg_exch_beta, "Matg_exch");
   }
   
   matx =NULL;
   memoria_double_uni(sizebi, &matx, "Matx");


// Creo arreglo para valores propios
 if (compara == 0) {
   valores =NULL;
   memoria_double_uni(sizedouble, &valores, "Valores");
 }
 else {
   valoresalfa =NULL;
   memoria_double_uni(sizedouble, &valoresalfa, "Valoresalfa");
   valoresbeta =NULL;
   memoria_double_uni(sizedouble, &valoresbeta, "Valoresbeta");
 }

// Creo matriz de KS depende el c\'alculo
//
    mat_ks =NULL;
    memoria_double_uni(sizebi, &mat_ks, "Mat_ks");
    for (i = 0; i < nt*nt; i++) mat_ks[i] = 0.f;
// Creo matriz de Fock

 if (compara == 0) {
   matfock =NULL;
   memoria_double_uni(sizebi, &matfock, "Matfock");
   for (i = 0; i < nt*nt; i++) matfock[i] = 0.f;
 }
 else {
   matfockalfa =NULL;
   memoria_double_uni(sizebi, &matfockalfa, "Matfockalfa");
   matfockbeta =NULL;
   memoria_double_uni(sizebi, &matfockbeta, "Matfockbeta");
   mat_ks_beta =NULL;
   memoria_double_uni(sizebi, &mat_ks_beta, "Mat_ks_beta");
   for (i = 0; i < nt*nt; i++) mat_ks_beta[i] = 0.f;
 }

// Creo arreglo de vectores finales

 if (compara == 0) {
     vectsfin =NULL;
     memoria_double_uni(sizebi, &vectsfin, "Vectoresfin");
 }
 else {
     vectsfinalfa =NULL;
     memoria_double_uni(sizebi, &vectsfinalfa, "Vectoresfinalfa");
     vectsfinbeta =NULL;
     memoria_double_uni(sizebi, &vectsfinbeta, "Vectoresfinbeta");
 }

// Creo arreglo de matriz de la densidad

 matp =NULL;
 memoria_double_uni(sizebi, &matp, "Matp");
 if (compara != 0) {
   matpalfa =NULL;
   memoria_double_uni(sizebi, &matpalfa, "Matpalfa");
   matpbeta =NULL;
   memoria_double_uni(sizebi, &matpbeta, "Matpbeta");
 }

// End of allocation memory

  printf("Memory to be allocated %8.2f Mb\n", (float) cuart*sizeof(double)/(1024*1024));

  for (i = 0; i < 100; i++) {
    arreglo_factorial[i] = factorial(i);
    arreglo_inv_factorial[i] = 1.f/arreglo_factorial[i];
  }

  if(strcmp(bound,"finite") == 0 || strcmp(bound,"dielectricc") == 0 || strcmp(bound,"polarization") == 0) {
     NC_minus = NULL;
     memoria_double_uni(sizedouble, &NC_minus, "NC_minus");
     NC_plus = NULL;
     memoria_double_uni(sizedouble, &NC_plus, "NC_plus");
     zetas = NULL;
     memoria_double_uni(sizedouble, &zetas, "zetas");
     if(strcmp(basis,"STOs") == 0 || strcmp(basis,"stos") == 0){
        for (i = 0; i < nt; i++){
          constants_normalization_finite(i, np, mang, expo, Rc,
                                         arreglo_factorial, arreglo_inv_factorial,
                                         using_gamma, gamma_couple, &tmp1, &tmp2, &tmp3, 0);
          NC_minus[i] = tmp1; //aquí se hace la asignación mike
          NC_plus[i] = tmp2;
          zetas[i] = tmp3;
        }
     }
  }

 char bound_pol[80];
 strcpy(bound_pol,"nothig");

 int iter_pol;
 double charge_int;

 charge_int = 0.f;
 iter_pol   = 0;
 
 potencial(using_gamma,
           nt, 
           z, 
           matv, 
           np, 
           mang, 
           ncm, 
           expo,
           basis,
           Rc, 
           gamma_couple,  
           bound, 
           iter_pol,
           charge_int,
           epsilon, 
           NC_minus, 
           NC_plus,
           arreglo_factorial, 
           arreglo_inv_factorial,
           elec,
           plasma);

 if (strcmp(bound,"polarization") == 0) {
    if(strcmp(basis,"STOs") == 0){
       sprintf(bound_pol,"%s", bound);
       strcpy(bound,"finite");
    }
 }

 if (strcmp(bound,"dielectricc") == 0) {
    if(strcmp(basis,"STOs") == 0){
      sprintf(bound_pol,"%s", bound);
      strcpy(bound,"finite");
    }
 }

 if (strcmp(bound,"dielectricnc") == 0) {
    if(strcmp(basis,"STOs") == 0){
      sprintf(bound_pol,"%s", bound);
      strcpy(bound,"free");
    }
 }


 traslape(using_gamma, 
          nt, 
          mats, 
          np, 
          mang, 
          ncm, 
          expo,
          basis, 
          Rc, 
          gamma_couple,
          bound, 
          epsilon, 
          NC_minus, 
          NC_plus, 
          arreglo_factorial,
          arreglo_inv_factorial); 

 cinetica(using_gamma,
          nt,  
          matk,  
          matkint,  
          matkext,  
          basis,
          np,  
          mang,  
          ncm,  
          expo,
          Rc,  
          gamma_couple,  
          bound,  
          NC_minus,  
          NC_plus,  
          arreglo_factorial,
          arreglo_inv_factorial);


 printf("Two-electron integrals on CPU\n");
 integrales_bie = (double *)malloc(cuart*sizeof(double *));
 if(integrales_bie == NULL) {
    fprintf(stderr, "out of memory\n");
    exit(1);
 }

 time_bie_ini = time(NULL);
 bielectronicas_CPU(using_gamma,
                          nt, 
                          np, 
                          mang, 
                          ncm, 
                          expo, 
                          bound, 
                          Rc, 
                          gamma_couple,
                          integrales_bie,
                          zetas, 
                          NC_minus, 
                          NC_plus, 
                          arreglo_factorial,
                          arreglo_inv_factorial);
 time_bie_fin = time(NULL);
 printf("Two-electron integrals (%g) in %.1f s\n", (float) cuart, (float) (time_bie_fin - time_bie_ini));
 
 int revisa, revisa1;
 revisa  = 0;
 revisa1 = 0;
 coef1 = (double)1;
 summat(dim, coef1, matk, coef1, matv, matcore);
 transforma(nt, mats, matx);
          // To start SCF
   if (compara == 0)
      copiamat(dim, matcore, matfock);
   else {
     copiamat(dim, matcore, matfockalfa);
     copiamat(dim, matcore, matfockbeta);
    }
   
   difmezcla = (double)1.0 - mezcla;
   energiavieja = (double)0.;
   iter = 0;     // mike, primera asignación de ITER
   iter_inter = 0;
   bandera = 0;
   coef1 = (double)1.;
   abre = 0;
   // zero = (double)0.;
   
   mat_temp_1 =NULL;
   memoria_double_uni(sizebi, &mat_temp_1, "Mat_temp_1");
   mat_temp_2 =NULL;
   memoria_double_uni(sizebi, &mat_temp_2, "Mat_temp_2");
   mat_temp_3 =NULL;
   memoria_double_uni(sizebi, &mat_temp_3, "Mat_temp_3");
   if (compara != 0) {
     mat_temp_1_beta =NULL;
     memoria_double_uni(sizebi, &mat_temp_1_beta, "Mat_temp_1_beta");
     mat_temp_2_beta =NULL;
     memoria_double_uni(sizebi, &mat_temp_2_beta, "Mat_temp_2_beta");
     mat_temp_3_beta =NULL;
     memoria_double_uni(sizebi, &mat_temp_3_beta, "Mat_temp_3_beta");
   }
   
   iter_inter = 0;
   bandera = 0;
   revisa = 0;
   revisa1 = 0;
   
   do {
       iter++;
       iter_pol++;
       if (compara == 0)  {
          for (i = 0; i < cuad; i++) {
             matg_exch[i] = (double) 0.0;
             matg_coul[i] = (double) 0.0;
             mat_ks[i]    = (double) 0.0;
          }

         if (strcmp(bound_pol,"polarization") == 0 && iter_pol > 1) {

            eigensystem_matdens(nt, tipo, elec, matx, matfock, valores,
                                vectsfin, matp, &revisa);

            grid_rhorad(z,
                        using_gamma,
                        compara,
                        bound,
                        basis,
                        nt,
                        elecalfa,
                        0,
                        Rc,
                        expo,
                        np,
                        zetas,
                        mang,
                        vectsfin,
                        NULL,
                        tipo,
                        NC_minus,
                        NC_plus,
                        gamma_couple,
                        grid,
                        grid_rho,
                        grid_rho_beta,
                        n_points,
                        save_i,
                        arreglo_factorial,
                        arreglo_inv_factorial);

               charge_int = 4.f*Pi*numerical_int(grid,
                                                 grid_rho,
                                                 save_i);

               printf("\n>>>>Elec  %3.4lf <<<<\n", charge_int);

               potencial(using_gamma,
                         nt,
                         z,
                         matv,
                         np,
                         mang,
                         ncm,
                         expo,
                         basis,
                         Rc,
                         gamma_couple,
                         bound_pol,
                         iter_pol,
                         charge_int,
                         epsilon,
                         NC_minus,
                         NC_plus,
                         arreglo_factorial,
                         arreglo_inv_factorial,
                         elec,
                         plasma);
    
              summat(dim, coef1, matk, coef1, matv, matcore);
              copiamat(dim, matcore, matfock);

              eigensystem_matdens(nt, tipo, elec, matx, matfock, valores,
                                  vectsfin, matp, &revisa);

           
         } else {
             eigensystem_matdens(nt, tipo, elec, matx, matfock, valores,
                                 vectsfin, matp, &revisa);
         }
       


  ///////COMIENZA_DFT
         if (flag_dft == 2 && strcmp(save_dft[1],"hf") == 0) {
           e_xc = 0.f;
         }
         else {
         grid_rhorad(z,
                     using_gamma,
                     compara,
                     bound,
                     basis,
                     nt,
                     elecalfa,
                     0,
                     Rc,
                     expo,
                     np,
                     zetas,
                     mang,
                     vectsfin,
                     NULL,
                     tipo,
                     NC_minus,
                     NC_plus,
                     gamma_couple,
                     grid,
                     grid_rho,
                     grid_rho_beta,
                     n_points,
                     save_i,
                     arreglo_factorial,
                     arreglo_inv_factorial);

         grid_derrad(z,
                     using_gamma,
                     compara,
                     bound,
                     basis,
                     nt,
                     elecalfa,
                     0,
                     Rc,
                     expo,
                     np,
                     zetas,
                     mang,
                     vectsfin,
                     NULL,
                     tipo,
                     NC_minus,
                     NC_plus,
                     gamma_couple,
                     grid,
                     grid_der,
                     grid_der_beta,
                     n_points,
                     save_i,
                     arreglo_factorial,
                     arreglo_inv_factorial);

         grid_secderrad(z,
                        using_gamma,
                        compara,
                        bound,
                        nt,
                        elecalfa,
                        0,
                        Rc,
                        expo,
                        np,
                        zetas,
                        mang,
                        vectsfin,
                        NULL,
                        tipo,
                        NC_minus,
                        NC_plus,
                        gamma_couple,
                        grid,
                        grid_secder,
                        grid_secder_beta,
                        n_points,
                        save_i,
                        arreglo_factorial,
                        arreglo_inv_factorial);

           e_xc = xc_energy(correlationc,
                            exchangex,
                            n_points,
                            compara,
                            save_dft,
                            weight_dft,
                            flag_dft,
                            array_i,
                            grid,
                            grid_rho,
                            NULL,
                            grid_der,
                            NULL,
                            grid_secder,
                            NULL);
      
          xc_potential(n_points,
                       compara,
                       save_i,
                       save_dft,
                       weight_dft,
                       flag_dft,
                       using_gamma,
                       bound,
                       Rc,
                       nt,
                       expo,
                       np,
                       zetas,
                       mang,
                       ncm,
                       NC_minus,
                       NC_plus,
                       gamma_couple,
                       array_i,
                       NULL,
                       array_ii,
                       NULL,
                       grid,
                       grid_rho,
                       NULL,
                       grid_der,
                       NULL,
                       grid_secder,
                       NULL,
                       mat_ks,
                       NULL,
                       arreglo_factorial,
                       arreglo_inv_factorial,
                       nolocal_pot_array);
         }

       } else { //ELSE COMPARA
           for (i = 0; i < cuad; i++) {
              matg_exch_alfa[i] = (double) 0.0;
              matg_exch_beta[i] = (double) 0.0;
           }

         if (strcmp(bound_pol,"polarization") == 0 && iter_pol > 1) {
           eigensystem_matdens(nt, tipo, elecalfa, matx, matfockalfa, valoresalfa,
                               vectsfinalfa, matpalfa, &revisa);
           eigensystem_matdens(nt, tipo, elecbeta, matx, matfockbeta, valoresbeta,
                               vectsfinbeta, matpbeta, &revisa1);
           coef1 = (double) 1.f;
           summat(dim, coef1, matpalfa, coef1, matpbeta, matp);

           grid_rhorad(z,
                       using_gamma,
                       compara,
                       bound,
                       basis,
                       nt,
                       elecalfa,
                       elecbeta,
                       Rc,
                       expo,
                       np,
                       zetas,
                       mang,
                       vectsfinalfa,
                       vectsfinbeta,
                       tipo,
                       NC_minus,
                       NC_plus,
                       gamma_couple,
                       grid,
                       grid_rho,
                       grid_rho_beta,
                       n_points,
                       save_i,
                       arreglo_factorial,
                       arreglo_inv_factorial);
               for(h = 0; h < n_points; h++) {
                  array_ii[h] = 0.f;
                  array_ii[h] = grid_rho[h] + grid_rho_beta[h];
               }



               charge_int = 4.f*Pi*numerical_int(grid,
                                                 array_ii,
                                                 save_i);

               printf("\n>>>>Elec  %3.4lf <<<<\n", charge_int);

               potencial(using_gamma,
                         nt,
                         z,
                         matv,
                         np,
                         mang,
                         ncm,
                         expo,
                         basis,
                         Rc,
                         gamma_couple,
                         bound_pol,
                         iter_pol,
                         charge_int,
                         epsilon,
                         NC_minus,
                         NC_plus,
                         arreglo_factorial,
                         arreglo_inv_factorial,
                         elec,
                         plasma);


             summat(dim, coef1, matk, coef1, matv, matcore);

             copiamat(dim, matcore, matfockalfa);
             copiamat(dim, matcore, matfockbeta);

             eigensystem_matdens(nt, tipo, elecalfa, matx, matfockalfa, valoresalfa,
                                 vectsfinalfa, matpalfa, &revisa);
             eigensystem_matdens(nt, tipo, elecbeta, matx, matfockbeta, valoresbeta,
                                 vectsfinbeta, matpbeta, &revisa1);
             coef1 = (double) 1.f;
             summat(dim, coef1, matpalfa, coef1, matpbeta, matp);


         } else {
            eigensystem_matdens(nt, tipo, elecalfa, matx, matfockalfa, valoresalfa,
                                vectsfinalfa, matpalfa, &revisa);
            eigensystem_matdens(nt, tipo, elecbeta, matx, matfockbeta, valoresbeta,
                                vectsfinbeta, matpbeta, &revisa1);
            coef1 = (double) 1.f;
            summat(dim, coef1, matpalfa, coef1, matpbeta, matp);
           }




           if (flag_dft == 2 && strcmp(save_dft[1],"hf") == 0) {
             e_xc = 0.f;
           }
           else {
               grid_rhorad(z,
                           using_gamma,
                           compara,
                           bound,
                           basis,
                           nt,
                           elecalfa,
                           elecbeta,
                           Rc,
                           expo,
                           np,
                           zetas,
                           mang,
                           vectsfinalfa,
                           vectsfinbeta,
                           tipo,
                           NC_minus,
                           NC_plus,
                           gamma_couple,
                           grid,
                           grid_rho,
                           grid_rho_beta,
                           n_points,
                           save_i,
                           arreglo_factorial,
                           arreglo_inv_factorial);

               grid_derrad(z,
                           using_gamma,
                           compara,
                           bound,
                           basis,
                           nt,
                           elecalfa,
                           elecbeta,
                           Rc,
                           expo,
                           np,
                           zetas,
                           mang,
                           vectsfinalfa,
                           vectsfinbeta,
                           tipo,
                           NC_minus,
                           NC_plus,
                           gamma_couple,
                           grid,
                           grid_der,
                           grid_der_beta,
                           n_points,
                           save_i,
                           arreglo_factorial,
                           arreglo_inv_factorial);

               grid_secderrad(z,
                              using_gamma,
                              compara,
                              bound,
                              nt,
                              elecalfa,
                              elecbeta,
                              Rc,
                              expo,
                              np,
                              zetas,
                              mang,
                              vectsfinalfa,
                              vectsfinbeta,
                              tipo,
                              NC_minus,
                              NC_plus,
                              gamma_couple,
                              grid,
                              grid_secder,
                              grid_secder_beta,
                              n_points,
                              save_i,
                              arreglo_factorial,
                              arreglo_inv_factorial);

                 e_xc = xc_energy(correlationc,
                                  exchangex,
                                  n_points,
                                  compara,
                                  save_dft,
                                  weight_dft,
                                  flag_dft,
                                  array_i,
                                  grid,
                                  grid_rho,
                                  grid_rho_beta,
                                  grid_der,
                                  grid_der_beta,
                                  grid_secder,
                                  grid_secder_beta);
            
                xc_potential(n_points,
                             compara,
                             save_i,
                             save_dft,
                             weight_dft,
                             flag_dft,
                             using_gamma,
                             bound,
                             Rc,
                             nt,
                             expo,
                             np,
                             zetas,
                             mang,
                             ncm,
                             NC_minus,
                             NC_plus,
                             gamma_couple,
                             array_i,
                             array_ii,
                             array_iii,
                             array_iv,
                             grid,
                             grid_rho,
                             grid_rho_beta,
                             grid_der,
                             grid_der_beta,
                             grid_secder,
                             grid_secder_beta,
                             mat_ks,
                             mat_ks_beta,
                             arreglo_factorial,
                             arreglo_inv_factorial,
                             nolocal_pot_array);
                }
         } //ELSE COMPARA
       if (revisa == 0 && revisa1 == 0) {
          //mrb  printf("\nSOY_YO2 %d %d\n", revisa, revisa1);
          matrixg_coul(nt, matp, matg_coul, np, mang, ncm, expo, Rc, bound,
                       integrales_bie);
          if (compara == 0) 
             matrixg_exch(nt, matp, matg_exch, np, mang, ncm, expo,
                          Rc, bound, tipo, integrales_bie);
          else {
              matrixg_exch(nt, matpalfa, matg_exch_alfa, np, mang, ncm, expo,
                           Rc, bound, tipo, integrales_bie);
              matrixg_exch(nt, matpbeta, matg_exch_beta, np, mang, ncm, expo,
                           Rc, bound, tipo, integrales_bie);
          }
            //Computing energy for each iteration
            //       compute_energy(compara, nt, matp, matk, matv, matg_coul,
            //                      &e_kin, &e_v, &e_core, &e_coul, &energia);
          if (iter > 7)
            energia_prueba = energia;

          compute_energy(compara, 
                         nt, 
                         matp, 
                         matk, 
                         matkint, 
                         matkext, 
                         matv, 
                         matg_coul, 
                         matg_exch,
                         matpalfa, 
                         matpbeta, 
                         matg_exch_alfa, 
                         matg_exch_beta,
                         &e_kin, 
                         &e_kinint, 
                         &e_kinext, 
                         &e_v, 
                         &e_core, 
                         &e_coul, 
                         &e_exch, 
                         &energia);
          //       energia = energia + e_exch;
         
        energia = energia  + coef_HF*e_exch + e_xc;

        if (flag_dft != 2) {
           if(compara == 0) {
             for (i = 0; i < cuad; i++) {
                temp_xcmat = 0.f;
                temp_xcmat = coef_HF*matg_exch[i] + mat_ks[i];
                matg_exch[i] = temp_xcmat;
             }
           } else {
               for (i = 0; i < cuad; i++) {
                 temp_xcmat = 0.f;
                 temp_xcmat = coef_HF*matg_exch_alfa[i] + mat_ks[i];
                 matg_exch_alfa[i] = temp_xcmat;
                 temp_xcmat = 0.f;
                 temp_xcmat = coef_HF*matg_exch_beta[i] + mat_ks_beta[i];
                 matg_exch_beta[i] = temp_xcmat;
               }
           }
        }

          if (iter <= 7)
            energia_prueba = energia;

          if (fabs(energia) < 2000.f && fabs(energia_prueba - energia) < 0.25) {
            printf("Iter=%5d---Energy=%14.10f", iter, energia);
            e_check = energia;
            difer = fabs(e_check - energiavieja);
          
            if (compara == 0) {
              //   DIIS para capa cerrada
              summat(dim, coef1, matg_coul, coef1, matg_exch, matg);
              summat(dim, coef1, matcore, coef1, matg, matfock);
              printf("\n");
              if (iter > 1000) {
                diis_check(nt, matp, mats, matfock, matx, mat_temp_2, &max_e);
                if ( max_e >= (double) 0.1) {
                    printf("\n");
                    iter_inter = 0;
                    bandera = 0;
                  } else {
                      printf("  DIIS = %10.6f\n", max_e);
                      bandera = 1;
                      iter_inter ++;
                      if (abre == 0) {
                       mat_f_store[iter_inter] = NULL;
                       memoria_double_uni(sizebi, &mat_f_store[iter_inter], "Mat_f_store");
                       mat_e_store[iter_inter] = NULL;
                       memoria_double_uni(sizebi, &mat_e_store[iter_inter], "Mat_e_store");
                      }
                      for ( i = 0; i < nt*nt; i++) {
                        mat_f_store[iter_inter][i] = matfock[i];
                        mat_e_store[iter_inter][i] = mat_temp_2[i];
                      }
                  }
              }
          
              if(iter > 1000) {
                if (bandera == 1 && iter_inter > 1) {
                  main_diis(nt, 
                              iter_inter, 
                              mat_temp_1, 
                              mat_temp_2, 
                              mat_temp_3,
                              mat_e_store, 
                              mat_f_store, 
                              matfock);
                } else {
                    summat(dim, difmezcla, matfock, mezcla, matgvieja, matfock);
                    copiamat(dim, matfock, matgvieja);
                  }
              }
                
            } else {// DIIS para capa abierta
                 summat(dim, coef1, matg_coul, coef1, matg_exch_alfa, matg_alfa);
                 summat(dim, coef1, matcore, coef1, matg_alfa, matfockalfa);
                 diis_check(nt, 
                            matpalfa, 
                            mats, 
                            matfockalfa, 
                            matx, 
                            mat_temp_2, 
                            &max_e);
                 summat(dim, coef1, matg_coul, coef1, matg_exch_beta, matg_beta);
                 summat(dim, coef1, matcore, coef1, matg_beta, matfockbeta);
                 diis_check(nt, 
                            matpbeta, 
                            mats, 
                            matfockbeta, 
                            matx, 
                            mat_temp_2_beta, 
                            &max_e_beta);
                 max_e_tope = (max_e > max_e_beta) ? max_e : max_e_beta;
                 if ( max_e_tope >= (double) 0.1) {
                   printf("\n");
                   iter_inter = 0;
                   bandera = 0;

                 } else {
                     printf("  DIIS = %10.6f, %10.6f\n", max_e, max_e_beta);
                     iter_inter ++;
                     if (abre == 0 ){
                     mat_f_store[iter_inter] = NULL;
                     memoria_double_uni(sizebi, &mat_f_store[iter_inter], "Mat_f_store");
                     mat_e_store[iter_inter] = NULL;
                     memoria_double_uni(sizebi, &mat_e_store[iter_inter], "Mat_e_store");
                     mat_f_store_beta[iter_inter] = NULL;
                     memoria_double_uni(sizebi, &mat_f_store_beta[iter_inter], "Mat_f_store");
                     mat_e_store_beta[iter_inter] = NULL;
                     memoria_double_uni(sizebi, &mat_e_store_beta[iter_inter], "Mat_e_store_beta");
                     }
                     for ( i = 0; i < nt*nt; i++) {
                       mat_f_store[iter_inter][i] = matfockalfa[i];
                       mat_e_store[iter_inter][i] = mat_temp_2[i];
                       mat_f_store_beta[iter_inter][i] = matfockbeta[i];
                       mat_e_store_beta[iter_inter][i] = mat_temp_2_beta[i];
                     }
                     bandera = 1;
                 }
                 
                 if (iter > 0) {
                   if (bandera == 1 && iter_inter > 1) {
                     main_diis(nt, 
                               iter_inter, 
                               mat_temp_1, 
                               mat_temp_2, 
                               mat_temp_3,
                               mat_e_store, 
                               mat_f_store, 
                               matfockalfa);
                     main_diis(nt, 
                               iter_inter, 
                               mat_temp_1_beta, 
                               mat_temp_2_beta, 
                               mat_temp_3_beta,
                               mat_e_store_beta, 
                               mat_f_store_beta, 
                               matfockbeta);
                   } else {
                       summat(dim, difmezcla, matg_alfa, mezcla, matgvieja_alfa, matg_alfa);
                       summat(dim, difmezcla, matg_beta, mezcla, matgvieja_beta, matg_beta);
                       copiamat(dim, matg_alfa, matgvieja_alfa);
                       copiamat(dim, matg_beta, matgvieja_beta);
                     }
                 }
                 if (iter_inter == 4) {
                    iter_inter = 0;
                    abre = 1;
                 }
              } //DIIS para capa abierta
                
            energiavieja = e_check;
           } else {
                printf("Iter=%5d---Energy=%14.8f. BlowUP!!\n", iter, energia);
                iter   = 1e7;
                difer  = tol;
                revisa = 1;
             }
         } else {
             iter   = 1e7;
             difer  = tol;
             revisa = 1;
           }
   } while(difer >= tol && iter <= maxiter);
   
// Writting XC potential. Check it!
   if (revisa == 0 && revisa1 == 0) {
//jgo     char name_out[50];
//jgo     if(flag_dft != 2) {
//jgo       sprintf(name_out,"Vxc_%s_%3.5lf_%3.5lf", bound, Rc, epsilon);
//jgo       write_out = fopen(name_out,"w");
//jgo       for (i = 0; i < n_points; i++)
//jgo         fprintf(write_out,"%16.12f %16.12f \n", grid[i], grid_rho[i]);
//jgo       fclose(write_out);
//jgo     }

             if(compara == 0) {
                double v_e;
                v_e = 0;
                for(i = 0; i < (elec/2); i++)
                   v_e = v_e + valores[i];
                printf("Var Energy                           = %15.5f\n", e_kin - e_v - 2.f*e_coul - 2.f*v_e);
             }

             if(strcmp(bound,"free") == 0 || strcmp(bound,"FREE") == 0){
                printf("Kinetic energy                       = %15.5f\n", e_kin);
                printf("Nuc-elec energy                      = %15.5f\n", e_v);
             }
             else{
                printf("Kinetic energy                       = %15.5f\n", e_kin);
                printf("Kinetic energy internal              = %15.5f\n", e_kinint);
                printf("Kinetic energy external              = %15.5f\n", e_kinext);     
                printf("Potential energy (V_int + V_ext)     = %15.5f\n", e_v);
             }
             printf("One-electron energy                  = %15.5f\n", e_core);
             printf("Coulomb energy                       = %15.5f\n", e_coul);
  //QUITAR DESPUES
             if (compara == 0) {
               if(strcmp(save_dft[1],"hf") == 0 && flag_dft == 2)
                 printf("Exchange energy                      = %15.5f\n", coef_HF*e_exch);
               else  {
                 printf("Exchange-Correlation energy          = %15.5f\n", coef_HF*e_exch + e_xc);
                 if (coef_HF*e_exch != 0)
                 printf("Exact-Exchange energy                = %15.5f\n", coef_HF*e_exch);
                 
                 for (i = 1; i < flag_dft; i++) { 
                 if (exchangex[i] != 0)               
                 printf("Functional %8s Exchange         = %15.5f\n", save_dft[i], exchangex[i]);
                 if (correlationc[i] != 0)
                 printf("Functional %8s Correlation      = %15.5f\n", save_dft[i], correlationc[i]);
                 }
                  
               }
             } 
             else
                printf("Exchange energy                      = %15.5f\n", e_exch);
             if(strcmp(save_dft[1],"hf") == 0 && flag_dft == 2) 
                printf("Hartree-Fock energy                  = %15.5f (%15.5f Ryd)\n", energia, 2.f*energia);
             else 
                printf("DFT energy                           = %15.5f (%15.5f Ryd)\n", energia, 2.f*energia);
 
            
             printf("POT/KIN                              = %15.5f\n", (e_v + e_coul + e_exch)/e_kin);
           
             if (iter >= maxiter) {
              printf("No convergence!!!!\n");
              valor = 1;     // mike
             }
             *total_energy = energia;     // mike
             time_scf_fin = time (NULL);
             printf("SCF time  = %16ld s.\n", time_scf_fin - time_scf_ini);
           
             printf("------------------- EIGENVALUES -------------------\n");
             if( compara == 0 ) {
                for(i = 0; i < nt; i++) 
                   printf("Eigenvalue %d: %8.5f\n", i, valores[i]);

                if(strcmp(properties,"property") == 0) {
                 printf("jgo, in properties\n");
                 xc_over_grid(compara, save_dft, flag_dft, weight_dft, n_points, grid,
                              grid_rho, grid_der);
                 wf_closed_shell(z, using_gamma, compara, bound, nt, elecalfa, elecbeta,
                                 Rc, expo, np, zetas, mang, ncm, vectsfin, NULL,
                                 tipo, NC_minus, NC_plus, gamma_couple,
                                 grid, grid_rho, NULL, grid_der, NULL, n_points,
                                 arreglo_factorial, arreglo_inv_factorial,
                                 matp, mats, &iter, save_i, print_vectors, cusp_kato);
                   print_out_array(n_points, grid, array_i, "xc.dat");
                   Evaluate_Elect_Pot(z, nt, matp, np, mang, ncm, expo, bound,
                                      arreglo_factorial, arreglo_inv_factorial, 
                                      grid, n_points, Rc, NC_minus, NC_plus, basis);
                }
             } 
             else { // Section for open-shell atoms
                for(i = 0; i < nt; i++) 
                   printf("Eigenvalue %d: alpha | beta = %8.5f  %8.5f\n", i, valoresalfa[i], valoresbeta[i]);
                grid_rhorad(z,
                            using_gamma,
                            compara,
                            bound,
                            basis,
                            nt,
                            elecalfa,
                            elecbeta,
                            Rc,
                            expo,
                            np,
                            zetas,
                            mang,
                            vectsfinalfa,
                            vectsfinbeta,
                            tipo,
                            NC_minus,
                            NC_plus,
                            gamma_couple,
                            grid,
                            grid_rho,
                            grid_rho_beta,
                            n_points,
                            save_i,
                            arreglo_factorial,
                            arreglo_inv_factorial);
               grid_derrad(z,
                           using_gamma,
                           compara,
                           bound,
                           basis,
                           nt,
                           elecalfa,
                           elecbeta,
                           Rc,
                           expo,
                           np,
                           zetas,
                           mang,
                           vectsfinalfa,
                           vectsfinbeta,
                           tipo,
                           NC_minus,
                           NC_plus,
                           gamma_couple,
                           grid,
                           grid_der,
                           grid_der_beta,
                           n_points,
                           save_i,
                           arreglo_factorial,
                           arreglo_inv_factorial);
            
           
                if (print_vectors == 1) {
                    int mu, nu, elemento1, elemento2;
                    double suma, rho_0, drho_0;
                    suma = 0.f;
                    for(mu = 0; mu < nt; mu++) {
                       for(nu = 0; nu < nt; nu++) {
                          elemento1 = mu*nt + nu;
                          elemento2 = nu*nt + mu;
                          suma = suma + matp[elemento1]*mats[elemento2];
                       }
                    }          
                    if(strcmp(bound,"free") == 0 || strcmp(bound,"finite") == 0 || strcmp(bound,"dielectricc") == 0 || strcmp(bound,"parabolic") == 0) {
//                    if(strcmp(bound,"finite") == 0) {     // esta era la sentencia original mike
                       printf("Number of electrons: = %f\n", suma);
                       printf("---------------------------\n");
                       printf("alpha\n");
                       expected_value_r(nt,
                                        -1,
                                        np,
                                        mang,
                                        ncm,
                                        expo,
                                        Rc,
                                        bound,
                                        basis,
                                        vectsfinalfa,
                                        elecalfa,
                                        gamma_couple,
                                        NC_minus,
                                        NC_plus,
                                        arreglo_factorial,
                                        arreglo_inv_factorial,
                                        using_gamma,
                                        zetas,
                                        grid);     
                       printf("beta\n");
                       expected_value_r(nt,
                                        -1,
                                        np,
                                        mang,
                                        ncm,
                                        expo,
                                        Rc,
                                        bound,
                                        basis,
                                        vectsfinbeta,
                                        elecbeta,
                                        gamma_couple,
                                        NC_minus,
                                        NC_plus,
                                        arreglo_factorial,
                                        arreglo_inv_factorial,
                                        using_gamma,
                                        zetas,
                                        grid);
                       printf("---------------------------\n");
                       printf("alpha\n");
                       expected_value_r(nt,
                                        1,
                                        np,
                                        mang,
                                        ncm,
                                        expo,
                                        Rc,
                                        bound,
                                        basis,
                                        vectsfinalfa,
                                        elecalfa,
                                        gamma_couple,
                                        NC_minus,
                                        NC_plus,
                                        arreglo_factorial,
                                        arreglo_inv_factorial,
                                        using_gamma,
                                        zetas,
                                        grid);
                       printf("beta\n");
                       expected_value_r(nt,
                                        1,
                                        np,
                                        mang,
                                        ncm,
                                        expo,
                                        Rc,
                                        bound,
                                        basis,
                                        vectsfinbeta,
                                        elecbeta,
                                        gamma_couple,
                                        NC_minus,
                                        NC_plus,
                                        arreglo_factorial,
                                        arreglo_inv_factorial,
                                        using_gamma,
                                        zetas,
                                        grid);
                       printf("---------------------------\n");
                       printf("alpha\n");
                       expected_value_r(nt,
                                        2,
                                        np,
                                        mang,
                                        ncm,
                                        expo,
                                        Rc,
                                        bound,
                                        basis,
                                        vectsfinalfa,
                                        elecalfa,
                                        gamma_couple,
                                        NC_minus,
                                        NC_plus,
                                        arreglo_factorial,
                                        arreglo_inv_factorial,
                                        using_gamma,
                                        zetas,
                                        grid);
                       printf("beta\n");
                       expected_value_r(nt,
                                        2,
                                        np,
                                        mang,
                                        ncm,
                                        expo,
                                        Rc,
                                        bound,
                                        basis,
                                        vectsfinbeta,
                                        elecbeta,
                                        gamma_couple,
                                        NC_minus,
                                        NC_plus,
                                        arreglo_factorial,
                                        arreglo_inv_factorial,
                                        using_gamma,
                                        zetas,
                                        grid);
                       printf("---------------------------\n");
                    }  //End if finite, free, dielec, parabolic

                   double SHAN;
                   for(h = 0; h < n_points; h++)
                      array_i[h] = 0.f;

                   array_i[0] = 0;
                   array_i[n_points - 1] = 0;
                   for(h = 1; h < n_points - 1; h++) {
                       if((grid_rho[h] + grid_rho_beta[h]) > 1e-20)
                        array_i[h] = -1.f*(grid_rho[h] + grid_rho_beta[h])*log((grid_rho[h] + grid_rho_beta[h]));
                   }

        
                   SHAN = 4.f*Pi*numerical_int(grid,
                                               array_i,
                                               n_points);

                   if(strcmp(bound,"free") == 0 && Rc == 0.f)
                     printf("\n%s %s Shannon entropy  %5.4lf\n", tipo, bound, SHAN);
                   else
                    printf("\n%s %s Shannon entropy %3.3lf %5.4lf\n", tipo, bound, Rc, SHAN);

                   for(h = 0; h < n_points; h++) {
                      array_ii[h] = 0.f;
                      array_ii[h] = grid_rho[h] + grid_rho_beta[h];
                   }

                   SHAN = 4.f*Pi*numerical_int(grid,
                                               array_ii,
                                               n_points);

                   printf("\n Number electrons from integral rho %5.8lf \n \n", SHAN);

                   FILE *workout;
                   char nameout[200];

                   if (strcmp(bound, "free") == 0 && Rc == 0.f) {
                    sprintf(nameout, "%s_%s_rho_drho_+drho_rdf_divdrhorho", bound, tipo);
                    if (strcmp(properties,"property") == 0) {
                      Evaluate_Elect_Pot(z, nt, matp, np, mang, ncm, expo, bound,
                                    arreglo_factorial, arreglo_inv_factorial,
                                    grid, n_points, Rc, NC_minus, NC_plus, basis);
                      print_out_array(n_points, grid, array_i, "xc.dat");
                    }
                   }
                   else
                    sprintf(nameout, "%3.3f_%s_%s_rho_drho_+drho_rdf_divdrhorho", Rc, bound, tipo);

                    workout = fopen(nameout, "w");

                   for (h = 0; h < n_points; h++)
                     if (grid_rho[h] + grid_rho_beta[h] > 1.e-16) 
                     fprintf(workout,"%5.4f  %24.14f  %24.14f  %24.14f  %24.14f %24.14f %24.14f\n", grid[h],
                                                                                                    grid_rho[h] + grid_rho_beta[h],
                                                                                                    grid_der[h] + grid_der_beta[h],
                                                                                                   -1.f*(grid_der[h] + grid_der_beta[h]),
                                                                                                    4.f*Pi*grid[h]*grid[h]*(grid_rho[h] + grid_rho_beta[h]),
                                                                                                    (grid_der[h] + grid_der_beta[h])/(grid_rho[h] + grid_rho_beta[h]),
                                                                                                   -1.f*(grid_der[h] + grid_der_beta[h])/(2.f*z*(grid_rho[h] + grid_rho_beta[h])));

                   fclose(workout);
          }//End for print_vectors
              
              
          printf("rho(0)                = %.4f\n", grid_rho[0] + grid_rho_beta[0]);
          printf("rho'(0)               = %.4f\n", grid_der[0] + grid_der_beta[0]);
          printf("KATO CUSP = %4.4f\n", -(grid_der[0] + grid_der_beta[0])/(2.f*z*(grid_rho[0] + grid_rho_beta[0])));
          *cusp_kato = -(grid_der[0] + grid_der_beta[0])/(2.f*z*(grid_rho[0] + grid_rho_beta[0]));
          if(strcmp(basis, "STOs") == 0 || strcmp(basis, "stos") == 0){
             if(-(grid_der[0] + grid_der_beta[0])/(2.f*z*(grid_rho[0] + grid_rho_beta[0]))  < 0.98 || -(grid_der[0] + grid_der_beta[0])/(2.f*z*(grid_rho[0] + grid_rho_beta[0])) > 1.20){
                iter = 1e7;
             //   printf("ojo \n"); mike: here was the problem
             }
          }
           ///////////////////PRUEBA DENSIDAD 
           
           
             fflush(stdout);
           //   }
            
           } // Finishing else for uhf 
           printf("---------------------------------------------------\n");
           
           
           if (iter >= maxiter) { 
             printf("Iter >= than maxiter, MP2 or EP2 cannot be done\n");
           } else { //label 1

                    mp2_energy = (double) 0.;
                    printf("\n%s %s\n", correlation, propagador);

                    if (strcmp(correlation,"mp2") == 0 || strcmp(propagador,"tom") == 0) { //label 3
                    //To consider important integrals 
                  //mrb  integrales = 0;
                  //mrb  tolbie = 1.e-7;
                  //mrb for (i = 0; i < cuart; i++)
                  //mrb   if (abs(integrales_bie[i]) > tolbie) {
                  //mrb     integrales++;
                  //mrb   }
           
                  //mrb   integral_importante = (int *)malloc(integrales*sizeof(int *));
                  //mrb//mrb if (strcmp(arq,"cpu") == 0)
                  //mrb   integrales_post_HF = (double *)malloc(integrales*sizeof(double *));
           
                  //mrb  j = -1;
                  //mrb for (i = 0; i < cuart; i++)
                  //mrb   if (abs(integrales_bie[i]) > tolbie) {
                  //mrb     j++;
                  //mrb     integral_importante[j] = i;
                  //mrb     integrales_post_HF[j] = integrales_bie[i];
                  //mrb   }
           
                  //mrb//mrb if (strcmp(arq,"cpu") == 0) {
                      vect_alfa_post = (double *)malloc(cuad*sizeof(double *));
                      vect_beta_post = (double *)malloc(cuad*sizeof(double *));
                   //mrb }
                
                    if (compara == 0) {
                      valoresalfa =NULL;
                      memoria_double_uni(sizedouble, &valoresalfa, "Valoresalfa");
                      valoresbeta =NULL;
                      memoria_double_uni(sizedouble, &valoresbeta, "Valoresbeta");
                      for (i = 0; i < nt; i++) {
                        valoresalfa[i] = valores[i];
                        valoresbeta[i] = valores[i];
                      }
                      for (i = 0; i < nt*nt; i++) {
                        vect_alfa_post[i] = vectsfin[i];
                        vect_beta_post[i] = vectsfin[i];
                      }
                    } else {
                        for (i = 0; i < nt*nt; i++) {
                          vect_alfa_post[i] = vectsfinalfa[i];
                          vect_beta_post[i] = vectsfinbeta[i];
                        }
                    }
           
                  //Starting MP2
                    if (strcmp(correlation,"mp2") == 0) { //label 2
           
                      printf("MP2 by using CPUs..\n");
                      arreglo_mo = (double *) malloc (cuart * sizeof (double));
                      if (arreglo_mo == NULL) {
                         puts ("\nError al reservar memoria para escalar");
                         exit (0);
                      }
           
                        mp2_cpu( nt,
                                 compara,
                                 elecalfa,
                                 elecbeta,
                                 vect_alfa_post,
                                 vect_beta_post,
                                 valoresalfa,
                                 valoresbeta,
                                 integrales_bie,
                                 arreglo_mo,
                                 &mp2_energy);
                   *total_energy += mp2_energy;
                    }//Finishing MP2 label 2
                   printf("Total energy           = %15.7f, %15.7f\n", *total_energy, energia);
           
                   if (strcmp(propagador,"tom") == 0) { //label 4
                      printf("EP2..\n");
                      arreglo_mo = (double *) malloc (cuart * sizeof (double));
                      if (arreglo_mo == NULL) {
                         puts ("\nError al reservar memoria para escalar");
                         exit (0);
                      }
                      arreglo_mo_ab = (double *) malloc (cuart * sizeof (double));
                      if (arreglo_mo_ab == NULL) {
                         puts ("\nError al reservar memoria para escalar");
                         exit (0);
                        }
                      ep2_cpu(espin,
                              nt,
                              elecalfa,
                              elecbeta,
                              integrales_bie,
                              vect_alfa_post,
                              vect_beta_post,
                              valoresalfa,
                              valoresbeta,
                              orbital,
                              arreglo_mo,
                              arreglo_mo_ab,
                              &omega);
                      *total_energy = omega;
                      free(arreglo_mo_ab);
                      free(arreglo_mo);
                   } //Finishing EP2 label 4
           
           
           
                   fflush(stdout);
                   if (compara == 0) {
                   free(valoresbeta);
                   free(valoresalfa);
                   } 
                  //mrb else {
                  //mrb free(integrales_post_HF);
                  //mrb }
                   free(vect_beta_post);
                   free(vect_alfa_post);
                  //mrb  free(integral_importante);
           
                   } //Finishing MP2 and EP2  label 3
           
           
            }//Finishing else maxiter label 1  
           
           
           
           // /////////////////////////////////////////////////////////////////////
           // Liberacion de memoria de todos los arreglos
           
          if (abre == 1) { 
           for (i = iter_inter; i >= 1; i--) {
             if (compara != 0) {
               free(mat_e_store_beta[i]);
               free(mat_f_store_beta[i]);
              } else { 
                  free(mat_e_store[i]);
                  free(mat_f_store[i]);
                }
           }
          } 
           
           if (compara != 0) {
             free(mat_temp_3_beta);
             free(mat_temp_2_beta);
             free(mat_temp_1_beta);
           }
           free(mat_temp_3);
           free(mat_temp_2);
           free(mat_temp_1);
           free(integrales_bie);
           free(matp);
           matp = 0;
           if (strcmp(bound,"finite") == 0) {
            free(zetas);
            free(NC_plus);
            free(NC_minus);
           }
           if (compara != 0) {
             free(matpalfa);
             matpalfa = 0;
             free(matpbeta);
             matpbeta = 0;
           }
           
           if (compara == 0 ) {
             free(vectsfin);
             vectsfin = 0;
           } else {
             free(vectsfinbeta);
             vectsfinbeta = 0;
             free(vectsfinalfa);
             vectsfinalfa = 0;
           }
           
           if (compara == 0) {
             free(matfock);
             matfock = 0;
             free(mat_ks);
             mat_ks = 0;
           } else {
             free(matfockbeta);
             matfockbeta = 0;
             free(matfockalfa);
             matfockalfa = 0;
           }
           
           if (compara == 0) {
             free(valores);
             valores = 0;
           } else {
             free(valoresbeta);
             valoresbeta = 0;
             free(valoresalfa);
             valoresalfa = 0;
           }
           
             free(matx);
             matx = 0;
             if (compara == 0) {
               free(matg_exch);
               matg_exch = 0;
               free(matg);
               matg = 0;
               free(matgvieja);
               matgvieja = 0;
             } else {
               free(matg_exch_beta);
               matg_exch_beta = 0;
               free(matg_exch_alfa);
               matg_exch_alfa = 0;
               free(matg_beta);
               matg_beta = 0;
               free(matg_alfa);
               matg_alfa = 0;
               free(matgvieja_beta);
               matgvieja_beta = 0;
               free(matgvieja_alfa);
               matgvieja_alfa = 0;
             }
             free(matg_coul);
             matg_coul = 0;
             free(matcore);
             matcore = 0;
             free(matv);
             matv = 0;
             free(matkext);
             matkext = 0;
             free(matkint);
             matkint = 0;
             free(matk);
             matk = 0;
             free(mats);
             mats = 0;
             free(nolocal_pot_array);
             nolocal_pot_array = 0;
             free(array_iv);
             array_iv = 0;
             free(array_iii);
             array_iii = 0;
             free(array_ii);
             array_ii = 0;
             free(array_i);
             array_i = 0;
           
             free(grid_secder_beta);
             grid_secder_beta = 0;
             free(grid_secder);
             grid_secder = 0;
           
             free(grid_der_beta);
             grid_der_beta = 0;
             free(grid_der);
             grid_der = 0;
             free(grid_rho_beta);
             grid_rho_beta = 0;
             free(grid_rho);
             grid_rho = 0;
             free(grid);
             grid = 0;

           
             time_4 = time (NULL);
             printf("Total time = %16ld s.\n", time_4 - time_3);
           
           /*mrb
           //Remover archivo vectores
            if ( remove(name_file_vectors) == 0 ) printf("\n has been removed\n");
            else printf("\n can't be removed\n");
           */
           
           //mrb vectsfin=NULL;
   } else {
           free(mat_temp_3);
           free(mat_temp_2);
           free(mat_temp_1);
           free(integrales_bie);
           free(matp);
           matp = 0;
           if (strcmp(bound,"finite") == 0) {
            free(zetas);
            free(NC_plus);
            free(NC_minus);
           }
           if (compara != 0) {
             free(matpalfa);
             matpalfa = 0;
             free(matpbeta);
             matpbeta = 0;
           }

           if (compara == 0 ) {
             free(vectsfin);
             vectsfin = 0;
           } else {
             free(vectsfinbeta);
             vectsfinbeta = 0;
             free(vectsfinalfa);
             vectsfinalfa = 0;
           }

           if (compara == 0) {
             free(valores);
             valores = 0;
           } else {
             free(valoresbeta);
             valoresbeta = 0;
             free(valoresalfa);
             valoresalfa = 0;
           }

             free(matx);
             matx = 0;
             if (compara == 0) {
               free(matg_exch);
               matg_exch = 0;
               free(matg);
               matg = 0;
               free(matgvieja);
               matgvieja = 0;
             } else {
               free(matg_exch_beta);
               matg_exch_beta = 0;
               free(matg_exch_alfa);
               matg_exch_alfa = 0;
               free(matg_beta);
               matg_beta = 0;
               free(matg_alfa);
               matg_alfa = 0;
               free(matgvieja_beta);
               matgvieja_beta = 0;
               free(matgvieja_alfa);
               matgvieja_alfa = 0;
             }
             free(matg_coul);
             matg_coul = 0;
             free(matcore);
             matcore = 0;
             free(matv);
             matv = 0;
             free(matkext);
             matkext = 0;
             free(matkint);
             matkint = 0;
             free(matk);
             matk = 0;
             free(mats);
             mats = 0;
             free(array_iv);
             array_iv = 0;
             free(array_iii);
             array_iii = 0;
             free(array_ii);
             array_ii = 0;
             free(array_i);
             array_i = 0;
           
             free(grid_secder_beta);
             grid_secder_beta = 0;
             free(grid_secder);
             grid_secder = 0;
           
             free(grid_der_beta);
             grid_der_beta = 0;
             free(grid_der);
             grid_der = 0;
             free(grid_rho_beta);
             grid_rho_beta = 0;
             free(grid_rho);
             grid_rho = 0;
             free(grid);
             grid = 0;
     }

    if(strcmp(bound_pol,"dielectricc") == 0 || strcmp(bound_pol,"polarization") == 0 || strcmp(bound_pol,"dielectricnc") == 0)
       sprintf(bound,"%s",bound_pol);

    if(iter >= maxiter || energia != energia){ // This is the original
       return(1);                                 //       se asigna que valga 1.
    } 
    else{  
       return(valor);
    }  
 }   // Termina scf

int exact_exch_ener(int     compara, 
                    int     nt, 
                    double *matp, 
                    double *matg_exch,
                    double *matpalfa, 
                    double *matpbeta,
                    double *matg_exch_alfa, 
                    double *matg_exch_beta,
                    double *e_exch)
{
 int i, j, element1, element2;
 double suma3, doble1;

 suma3 = (double)0.;

 if (compara == 0) {
   for (i = 0; i < nt; i++)
     for (j = 0; j < nt; j++) {
       element1 = j*nt + i;
       element2 = i*nt + j;
       suma3 += (double)0.5*matp[element2]*matg_exch[element1];
     }
  } else {
   for (i = 0; i < nt; i++)
     for (j = 0; j < nt; j++) {
       element1 = j*nt + i;
       element2 = i*nt + j;
       doble1 = matpalfa[element2]*matg_exch_alfa[element1];
       suma3 += (double)0.5*(doble1 + matpbeta[element2]*matg_exch_beta[element1]);
     }
 }

 *e_exch = suma3;

  return(0);
 }

