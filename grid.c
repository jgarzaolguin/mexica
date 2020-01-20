#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

 int building_grid(int           z,
                   char         *bound,
                   double        Rc,
                   double       *grid,
                   int           n_points,
                   int          *save_i)
 {//function

 int unsigned i;
 double point_r;
 double p1, pi;
 int save_temp, damping;

 pi = 4.f*atan(1.f);
 save_temp = 0;

 damping = 64;
 grid[0] = 0.f;

 for (i = 1; i < n_points; i++) { //for n_points
    point_r = -5.f + (double) (i - 1)/damping;                         // maya tipo Froese Fischer
    p1 = exp(point_r)/z;               // rj                      // p1 = [e^(point_r)]/z

    point_r = 1.f - cos(pi*((double) i/(2.f*(n_points+1))));
    p1 = 200.f*point_r;

//    if((p1 == Rc && strcmp(bound,"finite") == 0) || (strcmp(bound,"dielectricc") == 0) || (strcmp(bound,"polarization") == 0)){
    if((p1 == Rc && strcmp(bound,"free") == 0 ||  strcmp(bound,"finite") == 0) || (strcmp(bound,"dielectricc") == 0) || (strcmp(bound,"polarization") == 0)){
      point_r = -5.f + (double) (i)/damping;                         // maya tipo Froese Fischer
      p1 = exp(point_r)/z;    
     }
      
    grid[i] = p1;

    if (save_temp == 0) // Ya definiste save_temp = 0 va a entrar este cíclo
//      if (strcmp(bound,"polarization") == 0 || strcmp(bound,"finite") == 0 || strcmp(bound,"dielectricc") == 0 || strcmp(bound,"confined") == 0) { //if finite confined
      if (strcmp(bound,"free") == 0 ||  strcmp(bound,"polarization") == 0 || strcmp(bound,"finite") == 0 || strcmp(bound,"dielectricc") == 0 || strcmp(bound,"confined") == 0) { //if finite confined
        if (p1 > Rc) {   //este if sólo se cumple en cuanto p1 es mayor que R_{c} 
          if (i%2 != 0) {  //Acuerdate que el % es el operador de módulo o resto de una división (i/2) en este caso para i impar
             i = i + 1;    //reasignación
             grid[i] = Rc;
             grid[i - 1] = 0.5f*(grid[i] + grid[i - 2]);
           } else
               grid[i] = Rc; //más simple para i par
          *save_i    = i;
           save_temp = 1;   //Aquí reasignas el save_temp para que ya no se haga el if de save_temp
        }
      }

 } //for n_points
 return 0;
 }

