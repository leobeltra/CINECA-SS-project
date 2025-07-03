/*
################################  @ MAIN @   #####################################
#
# 
# 
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

    int Nbin_max = 100000;
    int Dbin=100;
    int Nbin =0;

    double **BINNED_CORRELATION = (double**) calloc(Nbin_max,sizeof(double*));
    double *C_MEAN = (double*) calloc(N,sizeof(double*));
    double *VAR = (double*) calloc(N,sizeof(double*));

    FILE *corr_bin;
    corr_bin = fopen("../../data_analysis/corr_bin_hermite2.txt", "wt");
    if( corr_bin ==NULL ) {
        perror("Error in opening file");
        exit(1);
    }
    FILE *corr_mean;
    corr_mean = fopen("../../data_analysis/corr_mean_hermite2.txt", "wt");
    if( corr_mean ==NULL ) {
        perror("Error in opening file");
        exit(1);
    }
    FILE *variance;
    variance = fopen("../../data_analysis/variance_hermite2.txt", "wt");
    if( variance ==NULL ) {
        perror("Error in opening file");
        exit(1);
    }
    
    /*
    #############################  @ TERMALIZZAZIONE @   #############################
    #
    # Termalizzazione a caldo della catena di Markov: vario il cammino della 
    # particella fino a termalizzazione della catena di markov, il cammino generato 
    # sarà accettato con probabilità alla Boltzmann exp(-dS), una volta che la catena 
    # ha raggiunto la condizione di ergodicità.
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
                    /*c = c + (4*xx[k]*xx[k]-2)*(4*xx[(k+t_phys)%N]*xx[(k+t_phys)%N]-2);*/
                    c = c + (xx[k]*xx[k]-1)*(xx[(k+t_phys)%N]*xx[(k+t_phys)%N]-1);
                }

                c=c/N;
                BINNED_CORRELATION[Nbin][t_phys] = BINNED_CORRELATION[Nbin][t_phys]+c; 
            }
            
        }
        
        Nbin=Nbin+1;
    }
    
    for (Nbin = 0; Nbin < Nbin_max; Nbin++)
    {
        for (t_phys = 0; t_phys < N ; t_phys++)
        {
            fprintf(corr_bin, "%lf\n", BINNED_CORRELATION[Nbin][t_phys]);
        } 
    }
    
    for (t_phys=0; t_phys<N; t_phys++)
    {
        double c=0;

        for (Nbin = 0; Nbin < Nbin_max; Nbin++)
        {
            BINNED_CORRELATION[Nbin][t_phys] = BINNED_CORRELATION[Nbin][t_phys]/(double)Dbin;
            c = c + BINNED_CORRELATION[Nbin][t_phys];
        }

        c=c/Nbin_max;
        C_MEAN[t_phys]=c;
        fprintf(corr_mean, "%lf\n", C_MEAN[t_phys]);
    }

    
    /*
    ##################   @ DEVIAZIONE STANDARD DEL CORRELATORE @   ###################
    #
    # Calcolo della deviazione standard del correlatore per ogni tempo fisico 
    # e salvataggio su file.
    #
    ##################################################################################
    */

    for (t_phys=0; t_phys<N; t_phys++)
    {
        for (Nbin = 0; Nbin < Nbin_max; Nbin++)
        {
            VAR[t_phys]=VAR[t_phys]+((BINNED_CORRELATION[Nbin][t_phys]-C_MEAN[t_phys])*(BINNED_CORRELATION[Nbin][t_phys]-C_MEAN[t_phys]));
        }
        VAR[t_phys]=VAR[t_phys]/(Nbin_max-1);
        fprintf(variance, "%lf\n", sqrt(VAR[t_phys]));
    }
    
    fclose(corr_bin);
    fclose(corr_mean);
    fclose(variance);
    
    return 0;    
}

