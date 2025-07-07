#define DH_ISING_2D_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"

double dH_ising_2D(double s_new, int i, int j, double** ss)
{

    double dH_B = 2*B_field*ss[i][j];

    double dH_S = 2*J*ss[i][j]*(ss[i][(j+1)%N]+ss[i][(j+N-1)%N]+ss[(i+1)%N][j]+ss[(i+N-1)%N][j]);

    double dH_ising=dH_B+dH_S;

    return dH_ising;

}
