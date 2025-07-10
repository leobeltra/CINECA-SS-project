#define H_ISING_2D_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"

double H_ising_2D(double* ss)
{

    double H_B=0;
    double H_S=0;
    int i,j;

    double H_S_rows=0;
    double H_S_cols=0;
    /*this must be at least in main */

    for ( i = 0; i < N; i++)
    {
        for (j = 0; j < N/2; j++)
        {
            H_B=H_B+ss[i*(N/2)+j];
        }
    }

    H_B = -B_field*H_B;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N/2; j++)
        {
            H_S_rows=H_S_rows+ss[i*(N/2)+j]*ss[i*(N/2)+((j+1)%(N/2))];
        }

    }

    for (j = 0; j < N/2; j++)
    {
        for (i = 0; i < N; i++)
        {
            H_S_cols=H_S_cols+ss[i*(N/2)+j]*ss[((i+1)%N)*(N/2) + j];
        }

    }

    H_S = -J*(H_S_rows+H_S_cols);

    double H_ising = H_S+H_B;

    return H_ising;

}

double H_issing_2D()
{

    double H_B=0;
    double H_S=0;
    int i,j;

    double H_S_rows=0;
    double H_S_cols=0;
    /*this must be at least in main */

    for ( i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            H_B=H_B+ss[i][j];
        }
    }

    H_B = -B_field*H_B;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            H_S_rows=H_S_rows+ss[i][j]*ss[i][(j+1)%N];
        }

    }

    for (j = 0; j < N; j++)
    {
        for (i = 0; i < N; i++)
        {
            H_S_cols=H_S_cols+ss[i][j]*ss[(i+1)%N][j];
        }

    }

    H_S = -J*(H_S_rows+H_S_cols);

    double H_ising = H_S+H_B;

    return H_ising;

}
