#define F_JACK_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../../include/global.h"

double f_jack( char fun_form, int n_obs, double data[]){

    void (*fun_ptr)(double) = &fun_form;
    double c[]=(*fun_ptr)(data[]);
    return c;
}
/*Per invocare la funzione in qualcos'altro

In ogni caso mi serve una stringa 
così io sul main scrivo esplicitamente la funzione che viene chiamata*/

void (*fun_ptr)(double) = &sin; /*al posto di sin deve finirci il contenuto della stringa*/