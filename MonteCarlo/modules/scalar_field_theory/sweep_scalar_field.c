#define SWEEP_SCALAR_FIELD_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "../../include/random.h"
#include "../../include/global.h"
#include "../../include/scalar_field_action.h"


double sweep_scalar_field(double r1[], double r2[]) {

    double dS;
    double accettanza = 0;
    int i, j, k;
    double phi_new;
    double phi_old;



    for(i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {

            for(k=0; k<N; k++){

                double S_old = scalar_field_action();

                phi_old = phi[i][j][k];
                phi_new=phi[i][j][k]+2*D*(r1[(N*i*j)+k]-0.5);

                phi[i][j][k]=phi_new;
                double S_new = scalar_field_action();

                dS = S_new-S_old;

                if(exp(-dS)>r2[(N*i*j)+k])
                {
                    phi[i][j][k]=phi_new;
                    accettanza = accettanza + 1;
                } 
                else{

                    phi[i][j][k]=phi_old;
                }

            }
            
        }
    }

    accettanza = accettanza/(N*N*N);
    return accettanza;
}