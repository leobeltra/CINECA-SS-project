#define AZIONE_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"



double azione_anharmonic_creutz()

{

double S=0;
int i;
/*azione in cui il termine quadratico compare con un meno*/
    for (i=0; i<N; i++){
        S  = S + ((M/2.0)*(xx[(i+1)%N]-xx[i])*(xx[(i+1)%N]-xx[i]) -A*(xx[i])*(xx[i]) + B*(xx[i]*xx[i]*xx[i]*xx[i]));
    }

return S;

}
