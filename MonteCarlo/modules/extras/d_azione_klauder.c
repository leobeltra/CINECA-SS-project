#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"
#include "../../include/random.h"
#include "../../include/azione.h"

double d_azione_klauder(double xxnew, int j){

    double delta_S ;
    double s = 0.817561;

    delta_S = (1./s)*(
                -0.5*M*(xx[(j+1)%N]-xxnew)*(xx[(j+1)%N]-xxnew)-0.5*M*(xxnew-xx[(j+N-1)%N])*(xxnew-xx[(j+N-1)%N])
    
                +0.5*M*(xx[(j+1)%N]-xx[j])*(xx[(j+1)%N]-xx[j])+0.5*M*(xx[j]-xx[(j+N-1)%N])*(xx[j]-xx[(j+N-1)%N])
              )

        +0.5*M*W*W*xx[j]*xx[j]-0.5*M*W*W*xxnew*xxnew

        - G*s*M*W*W*(xxnew*xxnew-xx[j]*xx[j]);


    /*delta_S = -0.5*M*(xx[(j+1)%N]-xxnew)*(xx[(j+1)%N]-xxnew)-0.5*M*(xxnew-xx[(j+N-1)%N])*(xxnew-xx[(j+N-1)%N])
    
                +0.5*M*(xx[(j+1)%N]-xx[j])*(xx[(j+1)%N]-xx[j])+0.5*M*(xx[j]-xx[(j+N-1)%N])*(xx[j]-xx[(j+N-1)%N])
                
                +0.5*M*M*G*((xx[(j+1)%N]-xx[j])*(xx[(j+1)%N]-xx[j])*(xx[(j+1)%N]-xx[j])*(xx[(j+1)%N]-xx[j]))

                +G*W*W;*/
                

        return delta_S;

}
