

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/azione.h"
#include "../../include/global.h"

int main() {


double s;
int j;

    for(j=0; j<N; j++){
        
        xx[j]=exp(j);
    }

s = azione();
printf("il valore dell'azione è %f", s );

return 0;

}
