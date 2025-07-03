#define CORRELATION_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"
#include "../../include/sweep.h"
#include "../../include/random.h"
#include "../../include/azione.h"

/*funzione che restituisce un double
prendendo in input t_markov e t_phys*/

double correlation ( int t_phys){
    int k;
    double c=0;
    for (k=0; k<N; k++)
    {
        c = c + xx[k]*xx[(k+t_phys)%N];
    }
    c=c/N;
    return c;
}
