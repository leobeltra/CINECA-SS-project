#define DAZIONE_AN_CREUTZ_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"
#include "../../include/random.h"
#include "../../include/azione.h"


double dazione_an_creutz(double xxnew, int j){

    double delta_S ;

    delta_S = +0.5*M*(xx[(j+1)%N]-xxnew)*(xx[(j+1)%N]-xxnew)+0.5*M*(xxnew-xx[(j+N-1)%N])*(xxnew-xx[(j+N-1)%N])
    
    -0.5*M*(xx[(j+1)%N]-xx[j])*(xx[(j+1)%N]-xx[j])-0.5*M*(xx[j]-xx[(j+N-1)%N])*(xx[j]-xx[(j+N-1)%N])

    -A*(xxnew*xxnew-xx[j]*xx[j]) 
    
    +B*((xxnew*xxnew*xxnew*xxnew)-(xx[j]*xx[j]*xx[j]*xx[j])) 
    
    /*+2*A*P*(xxnew-xx[j])*/

    ;

    return delta_S;


}
