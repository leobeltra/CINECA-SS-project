/*
############################  @ DECORRELATION @   ################################
#
# Federico De Matteis, Lab. di Fisica Computazionale 2022/2023
# Calcolo dell'autocorrelazione per la funzione di correlazione a due punti <x_lx_k>
# Binnaggio del correlatore e ricalcolo dell'autocorrelazione.
# 
##################################################################################
*/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/global.h"
#include "../include/sweep.h"
#include "../include/random.h"
#include "../include/dazione.h"
#include "../include/azione.h"

int main()
{
    int t_phys;
    int t_markov;
    double r3[N];
    double r4[N];

    int Nbin_max = 30000;
    int Dbin=40;
    int Nbin =0;
    int Dbin_max = 70;

    double **BINNED_CORRELATION = (double**) calloc(Nbin_max,sizeof(double*));
    double *GAMMA = (double*) calloc(N,sizeof(double*));


    FILE *deco;
    deco = fopen("../../data_analysis/deco.txt", "wt");
    if( deco ==NULL ) {
        perror("Error in opening file");
        exit(1);
    }
    /*
    #############################  @ TERMALIZZAZIONE @   #############################
    #
    # Termalizzazione a caldo della catena di Markov: vario il cammino della 
    # particella fino a termalizzazione della catena di markov, il cammino generato 
    # sarà distribuito con probabilità alla Boltzmann exp(-dS).
    #
    ##################################################################################
    */

    rlxd_init(2,251097);
    ranlxd(xx,N);

    for ( t_markov = 0; t_markov < 500; t_markov++)
    {
        ranlxd(r3,N);
        ranlxd(r4,N);
        sweep(r3,r4);
    }

    /*
    ##############   @ CORRELATORE BINNATO & MEDIA DEL CORRELATORE @   ###############
    #
    # Calcolo della matrice di correlazione binnata generando un vettore di 
    # correlazioni di N entrate per ogni sweep "; salvataggio della matrice del 
    # correlatore su file.
    #
    ##################################################################################
    */

    while (Nbin*Dbin < M_sweep){

        BINNED_CORRELATION[Nbin] = (double*) calloc(N,sizeof(double));
        for ( t_markov=Nbin*Dbin ; t_markov < (Nbin+1)*Dbin; t_markov++) 
        {
            ranlxd(r3,N);
            ranlxd(r4,N);
            sweep(r3,r4);

            for (t_phys = 0; t_phys < N ; t_phys++)
            {   
                int k;
                double c=0;
                
                for (k=0; k<N; k++)
                {
                    c = c + xx[k]*xx[(k+t_phys)%N];
                }

                c=c/N;
                BINNED_CORRELATION[Nbin][t_phys] = BINNED_CORRELATION[Nbin][t_phys]+c; 
                /*fprintf(corr_bin, "%lf\n", BINNED_CORRELATION[Nbin][t_phys]);*/
            }  
        }

        Nbin=Nbin+1;
    }

    /*
    ##########################   @ TEST DECORRELAZIONE  @   ##########################
    #
    # verifico che la covarianza del correlatore corrisponda a quella di variabili 
    # decorrelate e quindi indipendenti, calcolando il fattore Gamma(t_M)/Gamma(0)
    # per ogni tempo fisico e salvando i valori su file.
    # 
    ##################################################################################
    */

    for (t_phys =0; t_phys < 1; t_phys++) {

        double y_mean = 0;
        double Dbin_max = 70;

        for (t_markov = 0; t_markov < Nbin_max; t_markov++) {

            y_mean = y_mean + BINNED_CORRELATION[t_markov][t_phys];
        }
        y_mean= y_mean/(Nbin_max); 
    
        for (t_markov = 0; t_markov < Dbin_max; t_markov++) {

            double g = 0;
            double g_0 =0;
            int sweep;

            for (sweep = 0; sweep < Nbin_max-Dbin_max; sweep++) {

                g  = g + (BINNED_CORRELATION[sweep][t_phys]*BINNED_CORRELATION[sweep+t_markov][t_phys]);
                g_0 = g_0 +(BINNED_CORRELATION[sweep][t_phys]*BINNED_CORRELATION[sweep][t_phys]);
            }
        
            g = g/(Nbin_max-Dbin_max)-(y_mean*y_mean);
            g_0 = g_0/(Nbin_max-Dbin_max)-(y_mean*y_mean);
            
            fprintf(deco, "%lf\n", g/g_0);
    
        }
    }
    
    fclose(deco);
    return 0;
}
