#define DAZIONE_2D_AN_C
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"
#include "../../include/random.h"
#include "../../include/azione.h"

double dazione_2d_an(double xxnew, double yynew, int j){

    double delta_S ;


        delta_S = +0.5*M*(xx[(j+1)%N]-xxnew)*(xx[(j+1)%N]-xxnew)+0.5*M*(xxnew-xx[(j+N-1)%N])*(xxnew-xx[(j+N-1)%N])
    
        -0.5*M*(xx[(j+1)%N]-xx[j])*(xx[(j+1)%N]-xx[j])-0.5*M*(xx[j]-xx[(j+N-1)%N])*(xx[j]-xx[(j+N-1)%N])

        /*variazione dell'azione per il secondo gdl*/

        +0.5*M*(yy[(j+1)%N]-yynew)*(yy[(j+1)%N]-yynew)+0.5*M*(yynew-yy[(j+N-1)%N])*(yynew-yy[(j+N-1)%N])
    
        -0.5*M*(yy[(j+1)%N]-yy[j])*(yy[(j+1)%N]-yy[j])-0.5*M*(yy[j]-yy[(j+N-1)%N])*(yy[j]-yy[(j+N-1)%N])

        /*termine armonico*/

        -A*(xxnew*xxnew-xx[j]*xx[j])

        -A*(yynew*yynew-yy[j]*yy[j])

        /*termine anarmonico*/

        +B*(xxnew*xxnew*xxnew*xxnew-xx[j]*xx[j]*xx[j]*xx[j])

        +B*(yynew*yynew*yynew*yynew-yy[j]*yy[j]*yy[j]*yy[j])
        
        +2*B*(xxnew*xxnew*yynew*yynew-xx[j]*xx[j]*yy[j]*yy[j])
        ;


        return delta_S;

}
