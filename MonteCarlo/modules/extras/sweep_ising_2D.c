#define SWEEP_ISING_2D_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../../include/random.h"
#include "../../include/global.h"
#include "../../include/dH_ising_2D.h"

double sweep_ising_2D(double r[], double** ss) {

    double s_new;
    double accettanza = 0;
    int i, j;


    for(i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            s_new = -ss[i][j];

            if(exp(-beta*dH_ising_2D(s_new,i,j, ss))>r[(N)*i+j])
            {

                ss[i][j]=s_new;
                accettanza = accettanza + 1;
            } 
        }
    }

    accettanza = accettanza/(N*N);
    return accettanza;
}
