/* Integrales monoelectronicas y bielectronicas
   para resolver las ecuaciones de Hartree-Fock
   en sistemas atomicos libres.
   Jorge Garza, Junio del 2006 */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


double sto(int mu, double r, double* expo, int* np)
{
 double resconst, total;
 extern int constant_normalization_free(int* mu, double* expo, int* np, double* resultado);

 constant_normalization_free(&mu, expo, np, &resconst);


 if (r == 0) {
  if (np[mu] == 1)
    total = resconst;
  else
    total = 0.f;
 } else {
      total = resconst*exp(-expo[mu]*r);
     if (np[mu] == 1)
      total = total;
     else
      total = total*pow(r,(double) (np[mu] - 1));
     
   }

 return (total);
 }

double der_sto_r(int mu, double r, double* expo, int* np)
{
 double resconst, total;
 extern int constant_normalization_free(int* mu, double* expo, int* np, double* resultado);

 constant_normalization_free(&mu, expo, np, &resconst);

 if (r == 0.f) {
   if (np[mu] == 1)
     total = -resconst*expo[mu];
   else
   if (np[mu] == 2)
     total = resconst;
   else
     total = 0.f;
 } else {
     total = resconst*exp(-expo[mu]*r);

     if (np[mu] == 1)
       total = -total*expo[mu];
     else {
       total = total*pow(r,(double) (np[mu] - 1));
       total = total*((double) (np[mu] - 1)/r - expo[mu]);
     }
 }
 return (total);
 }

double sec_der_sto_r(int mu, double r, double* expo, int* np)
{
 double resconst, total;
 extern int constant_normalization_free(int* mu, double* expo, int* np, double* resultado);

 constant_normalization_free(&mu, expo, np, &resconst);

 if (r == 0.f) {
   if (np[mu] == 1)
     total = resconst*pow(expo[mu],2.f);
   else
   if (np[mu] == 2)
     total = -2.f*resconst*expo[mu];
   else
   if (np[mu] == 3)
     total = 2.f*resconst;
   else
     total = 0.f;
 } else {
     total = resconst*exp(-expo[mu]*r);

     if (np[mu] == 1)
       total = total*pow(expo[mu],2.f);
     else
     if (np[mu] == 2)
       total = total*(-2.f*expo[mu] + r*pow(expo[mu],2.f));
     else
     if (np[mu] == 3)
       total = total*(pow(r,2.f)*pow(expo[mu],2.f) + 2.f - 4.f*r*expo[mu]);
     else 
       total = total*(pow(r,(double) (np[mu] - 1))*pow(expo[mu],2.f) 
               + (double) pow(r,(double) (np[mu] - 3))*(np[mu] - 1)*(-2.f + np[mu] - 2.f*r*expo[mu]));
     
 }
 return (total);
 }


double free_orbital(int nt, int orbital, double r, double* expo, int* np,
                    double* vectors)
{
 int i;
  double suma;
  extern double sto(int mu, double r, double* expo, int* np);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*sto(i, r, expo, np);
  }

  return (suma);
 }

double der_free_orbital(int nt, int orbital, double r, double* expo, int* np,
                    double* vectors)
{
 int i;
  double suma;
  extern double der_sto_r(int mu, double r, double* expo, int* np);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*der_sto_r(i, r, expo, np);
  }

  return (suma);
 }

double sec_der_free_orbital(int nt, int orbital, double r, double* expo, int* np,
                    double* vectors)
{
 int i;
  double suma;
  extern double sec_der_sto_r(int mu, double r, double* expo, int* np);

  suma = 0.f;
  for (i = 0; i < nt; i++) {
    suma = suma + vectors[i*nt + orbital]*sec_der_sto_r(i, r, expo, np);
  }

  return (suma);
 }


double rho_radial_free(int nt, int elec, double r, double* expo,
                         int* np, double* vectors, char* tipo)
 {
  int orbital;
  double sum, partial, occ, pi;
  extern double free_orbital(int nt, int orbital, double r, double* expo, int* np,
                    double* vectors);

  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {
    partial = free_orbital(nt, orbital, r, expo, np, vectors);
    sum = sum + occ*partial*partial;
  }

  return(sum/(4.f*pi));
 }

double der_rho_radial_free(int nt, int elec, double r, double* expo,
                         int* np, double* vectors, char* tipo)
 {
  int orbital;
  double sum, partial, occ, pi;
  extern double free_orbital(int nt, int orbital, double r, double* expo, int* np,
                    double* vectors);
  extern double der_free_orbital(int nt, int orbital, double r, double* expo, int* np,
                    double* vectors);

  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {
    partial = 2.f*free_orbital(nt, orbital, r, expo, np, vectors);
    partial = partial*der_free_orbital(nt, orbital, r, expo, np, vectors);
    sum = sum + occ*partial;
  }

  return(sum/(4.f*pi));
 }

double sec_der_rho_radial_free(int nt, int elec, double r, double* expo,
                               int* np, double* vectors, char* tipo)
 {
  int orbital;
  double sum, partial, occ, pi;
  double partial2;
  extern double free_orbital(int nt, int orbital, double r, double* expo, int* np,
                    double* vectors);
  extern double der_free_orbital(int nt, int orbital, double r, double* expo, int* np,
                    double* vectors);
  extern double sec_der_free_orbital(int nt, int orbital, double r, double* expo, int* np,
                    double* vectors);

  pi = 4.f*atan(1.f);
  if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"RHF") == 0)
    occ = 2.f;
  else
    occ = 1.f;

  sum = 0.f;
  for (orbital = 0; orbital < elec; orbital++) {

    partial = 2.f*free_orbital(nt, orbital, r, expo, np, vectors);
    partial = partial*sec_der_free_orbital(nt, orbital, r, expo, np, vectors);

    partial2 = 2.f*der_free_orbital(nt, orbital, r, expo, np, vectors);
    partial2 = partial2*der_free_orbital(nt, orbital, r, expo, np, vectors);
    
    partial = partial + partial2;
   
    sum = sum + occ*partial;
  }

  return(sum/(4.f*pi));
 }

