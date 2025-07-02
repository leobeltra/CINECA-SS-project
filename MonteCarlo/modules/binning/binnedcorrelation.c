#define BINNEDCORRELATION_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"
#include "../../include/sweep.h"
#include "../../include/random.h"
#include "../../include/azione.h"


void binnedcorrelation(double **Y, double **C, int Dbin) {


    int t_markov;
    int t_phys;
    int Nbin;
    double c;

    /*questa funzione prende in input la matrice di correlazione e fa il binning
    è una funzione inutile perchè viene chiamata una sola volta*/

    for(t_phys=0; t_phys<N; t_phys ++)
    {
    
        Nbin=0;
        c=0;

        while(Dbin*Nbin < M_sweep){
        
            c=0;
            for (t_markov = Nbin*Dbin; t_markov<(Nbin+1)*Dbin; t_markov++){

                c = c + Y[t_phys][t_markov];
                
            }

            c = c/(double)Dbin;
            C[t_phys][Nbin] = c;
            
            Nbin=Nbin+1;
        }
        
    }

}

/*funzione che returna il puntatore alla testa del vettore di correlazioni*/

