#define DAZIONE_2D
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"
#include "../../include/random.h"
#include "../../include/azione.h"

double dazione_2d(double xxnew, double yynew, int j){

    double delta_S ;


        delta_S = -0.5*M*(xx[(j+1)%N]-xxnew)*(xx[(j+1)%N]-xxnew)-0.5*M*(xxnew-xx[(j+N-1)%N])*(xxnew-xx[(j+N-1)%N])
    
        -0.5*M*W*W*xxnew*xxnew+0.5*M*(xx[(j+1)%N]-xx[j])*(xx[(j+1)%N]-xx[j])+0.5*M*(xx[j]-xx[(j+N-1)%N])*(xx[j]-xx[(j+N-1)%N])

        +0.5*M*W*W*xx[j]*xx[j]

        /*variazione dell'azione per il secondo gdl*/

        -0.5*M*(yy[(j+1)%N]-yynew)*(yy[(j+1)%N]-yynew)-0.5*M*(yynew-yy[(j+N-1)%N])*(yynew-yy[(j+N-1)%N])
    
        -0.5*M*W*W*yynew*yynew+0.5*M*(yy[(j+1)%N]-yy[j])*(yy[(j+1)%N]-yy[j])+0.5*M*(yy[j]-yy[(j+N-1)%N])*(yy[j]-yy[(j+N-1)%N])

        +0.5*M*W*W*yy[j]*yy[j]

        +lambda*(xx[j]*xx[j]*xx[j]*xx[j]-xxnew*xxnew*xxnew*xxnew)

        +lambda*(yy[j]*yy[j]*yy[j]*yy[j]-yynew*yynew*yynew*yynew)
        ;


        return delta_S;

}
