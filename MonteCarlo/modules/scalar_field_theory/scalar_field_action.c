#define SCALAR_FIELD_ACTION_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "../../include/global.h"



double scalar_field_action() /*1.186/1.220 on Montvay-Munster*/

{
    int t, x, y;
    double S=0;
    double S_space=0;
    double S_time=0;
    double S_space_y=0;
    double _Complex i = csqrt(-1);
    double _Complex Z =  0.5 * i;
    /*printf("%f + i%f\n", creal(Z), cimag(Z));*/

    /*for(t=0; t<N; t++){
        
        S_time=S_time-2*K*phi[t][x][y]*phi[(t+1)%N][x][y]+phi[t][x][y]*phi[t][x][y]+G*(phi[t][x][y]*phi[t][x][y]-1)*(phi[t][x][y]*phi[t][x][y]-1)-G;
    }

    for(x=0; x<N; x++){

        S_space=S_space-2*K*phi[t][x][y]*phi[t][(x+1)%N][y]+phi[t][x][y]*phi[t][x][y]+G*(phi[t][x][y]*phi[t][x][y]-1)*(phi[t][x][y]*phi[t][x][y]-1)-G;
        
    }
    for(y=0; y<N; y++){

        S_space_y=S_space_y-2*K*phi[t][x][y]*phi[t][x][(y+1)%N]+phi[t][x][y]*phi[t][x][y]+G*(phi[t][x][y]*phi[t][x][y]-1)*(phi[t][x][y]*phi[t][x][y]-1)-G;
        
    }

        
    return S=S_space+S_time+S_space_y;
        
}*/
    


    /*time like part of the action-sum over time index*/
    for(t=0; t<N; t++){

        double L=0;
        double L_T=0;
        double L_T_1=0;

        for(x=0; x<N; x++){
            for(y=0; y<N; y++){

                L=L+0.5*(phi[(t+1)%N][x][y]-phi[t][x][y])*(phi[(t+1)%N][x][y]-phi[t][x][y]);

            }
        }

        for(x=0; x<N; x++){
            for(y=0; y<N; y++){
            L_T = L_T + 0.5*(phi[t][x][y]-phi[t][(x+1)%N][y])*(phi[t][x][y]-phi[t][(x+1)%N][y])
            +0.5*(phi[t][x][y]-phi[t][x][(y+1)%N])*(phi[t][x][y]-phi[t][x][(y+1)%N])
            +0.5*M*M*phi[t][x][y]*phi[t][x][y]+phi[t][x][y]*phi[t][x][y]*phi[t][x][y]*phi[t][x][y];
            }
        }

        for(x=0; x<N; x++){
            for(y=0; y<N; y++){
            L_T_1 = L_T_1 + 0.5*(phi[(t+1)%N][x][y]-phi[(t+1)%N][(x+1)%N][y])*(phi[(t+1)%N][x][y]-phi[(t+1)%N][(x+1)%N][y])
            +0.5*(phi[(t+1)%N][x][y]-phi[(t+1)%N][x][(y+1)%N])*(phi[(t+1)%N][x][y]-phi[(t+1)%N][x][(y+1)%N])
            +0.5*M*M*phi[(t+1)%N][x][y]*phi[(t+1)%N][x][y]+phi[(t+1)%N][x][y]*phi[(t+1)%N][x][y]*phi[(t+1)%N][x][y]*phi[(t+1)%N][x][y];
            }
        }

        S=S+L+0.5*L_T+0.5*L_T_1;
    }
    return S;
}

