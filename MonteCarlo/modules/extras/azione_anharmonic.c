#define AZIONE_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"



double azione_anharmonic()

{

double S=0;
int i;

    for (i=0; i<N; i++){
        S  = S + (M/2.0)*(xx[(i+1)%N]-xx[i])*(xx[(i+1)%N]-xx[i])+((M*W*W)/2.0)*(xx[i])*(xx[i]) + lambda*xx[i]*xx[i]*xx[i]*xx[i];
    }

return S;

}
