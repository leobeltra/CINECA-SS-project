#define H_ISING_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"

double H_Ising()

{

double H_B=0;
double H_S=0;
int i;

/*this must e at least in main */

for ( i = 0; i < N; i++)
{
    H_B=H_B+s[i];
}
H_B = -B_field*H_B;
for (i = 0; i < N; i++)
{
    H_S =H_S+s[i]*s[(i+1)%N];
}
H_S = -J*H_S;
double H_ising = H_S+H_B;

return H_ising;

}
