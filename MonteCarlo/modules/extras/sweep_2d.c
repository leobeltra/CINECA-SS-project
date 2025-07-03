#define SWEEP_2D_C

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
#include "../../include/dazione_2d.h"
#include "../../include/dazione_2d_an.h"

double sweep_2d(double r1[]) {

    double xxnew;
    double yynew;
    double accettanza = 0;
    int i;


    for(i=0; i<N; i++){
           
        xxnew = xx[i] - 2*D*(r1[i]-0.5);
        yynew = yy[i] - 2*E*(r1[i+N]-0.5);

        if(exp(-dazione_2d_an(xxnew,yynew,i))>r1[i+(2*N)]){

            xx[i]=xxnew;
            yy[i]=yynew;

            accettanza = accettanza + 1;
        } 
    }

    accettanza = accettanza/N;
    return accettanza;
}
