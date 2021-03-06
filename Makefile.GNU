# Autor: Jorge Garza
#
# To use several cores, you need the next instruction
#  export OMP_NUM_THREADS=16
# For this case, mexica will use 16 cores

CC = gcc

# In this makefile lapack and blas are installed in DIRLIB
DIRLIB=${HOME}/programs/lapack-3.5.0
librerias=-L$(DIRLIB) $(DIRLIB)/liblapack.a $(DIRLIB)/librefblas.a

objects = data.o optimiza.o clebsch.o integrals_free.o integrals_finite.o integrals_confined.o normalization_constants.o checkout_expotents_finite.o  main.o matrices.o traslape.o cinetica.o potencial.o scf.o matrixg_coul.o  matrixg_exch.o diis.o mp2.o basis_assigned.o basis_extern.o ep2.o shannon_entropy_finite.o shannon_entropy_confined.o main_function_global.o spherical_density_finite.o spherical_density_confined.o spherical_density_free.o print_spherical_density_finite.o check_linear_dependence.o

flags = -c -O1
mexica : $(objects)
	$(CC) -o mexica $(objects) $(librerias) -fopenmp -lm -L/usr/lib64 -lgfortran

data.o : data.c
	$(CC) $(flags) data.c
basis_assigned.o : basis_assigned.c
	$(CC) $(flags) basis_assigned.c
basis_extern.o : basis_extern.c
	$(CC) $(flags) basis_extern.c
optimiza.o : optimiza.c
	$(CC) $(flags) optimiza.c
main.o : main.c
	$(CC) $(flags) main.c
clebsch.o : clebsch.c
	$(CC) $(flags) clebsch.c
integrals_free.o : integrals_free.c
	$(CC) $(flags) integrals_free.c
integrals_finite.o : integrals_finite.c
	$(CC) $(flags) integrals_finite.c
integrals_confined.o : integrals_confined.c
	$(CC) $(flags) integrals_confined.c
normalization_constants.o : normalization_constants.c
	$(CC) $(flags) normalization_constants.c
checkout_expotents_finite.o : checkout_expotents_finite.c
	$(CC) $(flags) checkout_expotents_finite.c
matrices.o : matrices.c
	$(CC) $(flags) matrices.c
traslape.o : traslape.c
	$(CC) $(flags) -fopenmp -lpthread traslape.c
cinetica.o : cinetica.c
	$(CC) $(flags) -fopenmp -lpthread cinetica.c
potencial.o : potencial.c
	$(CC) $(flags) -fopenmp -lpthread potencial.c
eigensystem.o : eigensystem.c
	$(CC) $(flags) eigensystem.c
##nrutil.o : nrutil.c
##	$(CC) $(flags) nrutil.c
matrixg_coul.o : matrixg_coul.c
	$(CC) $(flags) -fopenmp -lpthread matrixg_coul.c
matrixg_exch.o : matrixg_exch.c
	$(CC) $(flags) -fopenmp -lpthread matrixg_exch.c
diis.o : diis.c
	$(CC) $(flags) diis.c
mp2.o : mp2.c
	$(CC) $(flags) -fopenmp -lpthread mp2.c
ep2.o : ep2.c
	$(CC) $(flags) ep2.c
scf.o : scf.c
	$(CC) $(flags) -fopenmp -lpthread scf.c
shannon_entropy_finite.o : shannon_entropy_finite.c
	$(CC) $(flags) shannon_entropy_finite.c 
shannon_entropy_confined.o : shannon_entropy_confined.c
	$(CC) $(flags) shannon_entropy_confined.c 
main_function_global.o : main_function_global.c
	$(CC) $(flags) main_function_global.c 
spherical_density_finite.o : spherical_density_finite.c
	$(CC) $(flags) spherical_density_finite.c
spherical_density_confined.o : spherical_density_confined.c
	$(CC) $(flags) spherical_density_confined.c
spherical_density_free.o : spherical_density_free.c
	$(CC) $(flags) spherical_density_free.c
print_spherical_density_finite.o : print_spherical_density_finite.c
	$(CC) $(flags) print_spherical_density_finite.c
check_linear_dependence.o : check_linear_dependence.c
	$(CC) $(flags) check_linear_dependence.c


##grid_ks.o : fortran/grid_ks.f
##	gfortran -c fortran/grid_ks.f
##gauss_quadrature.o : fortran/gauss_quadrature.f
##	gfortran -c fortran/gauss_quadrature.f
##matrixg_ks.o : fortran/matrixg_ks.f
##	gfortran -c fortran/matrixg_ks.f
##rho_derrho_on_mesh.o : rho_derrho_on_mesh.c
##	$(CC) $(flags) rho_derrho_on_mesh.c
clean : 
	rm $(objects) mexica
