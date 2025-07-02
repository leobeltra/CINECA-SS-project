#define AZIONE_BF_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"



double azione_bf()

{

double S=0;
int i;

    for (i=0; i<N; i++){
        S  = S
        
        + 0.5*(
                (M/2.0)*(xx[(i+1)%N]-xx[i])*(xx[(i+1)%N]-xx[i]) 
            +   (M/2.0)*(xx[(i)]-xx[(i+N-1)%N])*(xx[(i)]-xx[(i+N-1)%N] ) 
        )
        
        + ((M*W*W)/2.0)*(xx[i])*(xx[i]);
    }

return S;

}
