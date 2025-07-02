#define SWEEP_AN_CREUTZ_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../../include/random.h"
#include "../../include/global.h"
#include "../../include/sweep.h"
#include "../../include/dazione.h"
#include "../../include/dazione_an_creutz.h"

double sweep_an_creutz(double r1[],double r2[]) {

    double xxnew;
    double accettanza = 0;
    int i;


    for(i=0; i<N; i++){
           
        xxnew = xx[i] - 2*D*(r1[i]);

        if(exp(dazione_an_creutz(xxnew,i))>r2[i]){

            xx[i]=xxnew;
            accettanza = accettanza + 1;
        } 
    }

    accettanza = accettanza/N;
    return accettanza;
}

