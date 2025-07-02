#define GEOMETRIZE_V_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"


double geometrize_V(int start, double q, double xx_init, int j){

    int i =0;
    double power_q=1;

    for (i = start; i < 2*j; i++)
    {
        power_q = power_q*q;
    }

    double V_geom = (1/2.)*xx_init*xx_init*((1-power_q)/(1-q*q));

    return V_geom;
    
}