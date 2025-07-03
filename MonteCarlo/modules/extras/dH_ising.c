#define DH_ISING_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"

double dH_ising(double s_new, int j){

    double dH_ising;
    double dH_S = 2*J*s[j]*(s[(j+1)%N]+s[(j+N-1)%N]);
    double dH_B = 2*B_field*s[j];
    dH_ising = dH_B + dH_S;

    return dH_ising;

}