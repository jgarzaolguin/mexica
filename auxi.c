int grid_rhorad_orbital(int           z,
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
                int           selected_orb)
 
 {//function

 extern double rho_orbital_radial_finite_int(char   *using_gamma,
                                     int     nt,
                                     int     elec,
                                     double  r,
                                     double  Rc,
                                     double *zeta,
                                     int    *np,
                                     double *vectors,
                                     char   *tipo,
                                     double  gamma_couple,
                                     double *N_minus,
                                     int     selected_orb);

 extern double rho_orbital_radial_finite_ext(int     nt,
                                     int     elec,
                                     double  r,
                                     double  Rc,
                                     double *alfa,
                                     int    *mang,
                                     double *vectors,
                                     char   *tipo,
                                     double  gamma_couple,
                                     double *N_plus,
                                     int     selected_orb);

 int unsigned i;
 double p1;


 for (i = 0; i < n_boundary; i++) { //for save_temp
   p1 = grid[i];
   if (compara == 0) {
     grid_rho[i] = rho_orbital_radial_finite_int(using_gamma, nt, elecalfa, p1, Rc, expo, np,
                                    vectsfinalfa, tipo, gamma_couple, NC_minus, selected_orb);
   } else {
         grid_rho[i] = rho_orbital_radial_finite_int(using_gamma, nt, elecalfa, p1, Rc, expo, np,
                                        vectsfinalfa, tipo, gamma_couple, NC_minus, selected_orb);
         grid_rho_beta[i] = rho_orbital_radial_finite_int(using_gamma, nt, elecbeta, p1, Rc, expo, np,
                                             vectsfinbeta, tipo, gamma_couple, NC_minus, selected_orb);
   }
 } //for save_temp
         
 for (i = n_boundary; i < n_points; i++) { //for save_temp
   p1 = grid[i];
   if (compara == 0) {
     grid_rho[i] = rho_orbital_radial_finite_ext(nt, elecalfa, p1, Rc, zetas, mang, 
                           vectsfinalfa, tipo, gamma_couple, NC_plus, selected_orb);
   } else {
         grid_rho[i] = rho_orbital_radial_finite_ext(nt, elecalfa, p1, Rc, zetas, mang,
                               vectsfinalfa, tipo, gamma_couple, NC_plus, selected_orb);
         
         grid_rho_beta[i] =  rho_orbital_radial_finite_ext(nt, elecbeta, p1, Rc, zetas, mang,
                                     vectsfinbeta, tipo, gamma_couple, NC_plus, selected_orb);
     }
 }
 return 0;
}//function

