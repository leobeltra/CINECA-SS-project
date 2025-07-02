#define SWEEP_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../../include/random.h"
#include "../../include/global.h"
#include "../../include/sweep.h"
#include "../../include/dazione.h"
#include "../../include/dazione_an.h"
#include "../../include/dazione_an_creutz.h"
#include "../../include/d_azione_bf.h"
#include "../../include/d_azione_klauder.h"



double sweep(double r1[],double r2[]) {

    double xxnew;
    double accettanza = 0;
    int i;


    for(i=0; i<N; i++){
           
        xxnew = xx[i] - 2*D*(r1[i]-0.5);

        if(exp(dazione(xxnew,i))>r2[i]){

            xx[i]=xxnew;
            accettanza = accettanza + 1;
            /*printf(""exp(-dazione_an_creutz(xxnew,i)))*/
        } 
    }

    accettanza = accettanza/N;
    return accettanza;
}
