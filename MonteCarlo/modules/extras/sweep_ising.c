#define SWEEP_ISING_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../../include/random.h"
#include "../../include/global.h"
#include "../../include/dH_ising.h"

double sweep_ising(double r[], double** s) {

    double s_new;
    double accettanza = 0;
    int i;


    for(i=0; i<N; i++){
           
        s_new = -s[i];

        if(exp(-beta*dH_ising(s_new,i, s))>r[i]){

            s[i]=s_new;
            accettanza = accettanza + 1;
        } 
    }

    accettanza = accettanza/N;
    return accettanza;
}
