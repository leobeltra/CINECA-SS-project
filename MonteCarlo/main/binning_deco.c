#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/global.h"
#include "../include/sweep.h"
#include "../include/random.h"
#include "../include/azione.h"
#include "../include/binnedcorrelation.h"

int main() {
    /*initializing variables*/
    int t_real;
    int t_markov;
    int j;
    double y;
    int i;
    int k;
    double r3[N];
    double r4[N];
    int level=1;
    int seed=123678;
    
    int t_phys;
    double c;
    int Nbin=0;
    double Nbin_max = 30000;
    int Dbin=40;
    double c_bar=0;
    /*malloc matrix*/

    double **C = (double**) malloc(N*sizeof(double*)); 

    for(i=0; i<N; i++)
    {
	    C[i] = (double*) malloc(Nbin_max*sizeof(double));
    }

    double **Y = (double**) malloc(N*sizeof(double*)); 

    for(i=0; i<N; i++)
    {
	    Y[i] = (double*) malloc(M_sweep*sizeof(double));
    }

    /*Termalization of the markov chain*/
    /*cold initialization*/
    for ( j = 0; j<N; j++) {
        xx[j]=0;
    }

    void rlxd_init(int level,int seed);

    /*termalization of the markov chain using 10 times t_markov of termalization */
    for ( j = 0; j < 500; j++)
    {
        ranlxd(r3,N);
        ranlxd(r4,N);
        sweep(r3,r4);
    }

    FILE *file_C_binned;
    file_C_binned = fopen("../../data_analysis/C_binned.txt", "wt");
    if( file_C_binned ==NULL ) {
        perror("Error in opening file");
        exit(1);
    }

    FILE *file_C_ensemble;
    file_C_ensemble = fopen("../../data_analysis/C_mean_ensemble.txt", "wt");
    if( file_C_ensemble ==NULL ) {
        perror("Error in opening file");
        exit(1);
    }

    for ( t_markov = 0; t_markov < M_sweep; t_markov++) {

        ranlxd(r3,N);
        ranlxd(r4,N);
        sweep(r3,r4);

        /*compute Y[N]*/

        for (t_real = 0; t_real < N ; t_real++) {
        
            y=0;

            for (k=0; k<N; k++){
                y = y + xx[k]*xx[(k+t_real)%N];
            }

            y=y/N;
            Y[t_real][t_markov] = y; /*correlation matrix*/
        }
    }
    
    /*Binnaggio della funzione di correlazione
    non serve creare un fnzione per usarla una sola volta*/
    binnedcorrelation(Y,C,Dbin);

    for (t_phys=0; t_phys<N; t_phys ++)
    {
        free(Y[t_phys]);
    }
    

    double C_bar[N];
    /*averaging the correlation function in the ensemble*/
    for(t_phys=0; t_phys<N; t_phys++){

        for(Nbin=0; Nbin < Nbin_max; Nbin++){
            c_bar=c_bar + C[t_phys][Nbin];
            /*printf("%lf\n",C[t_phys][Nbin]);*/
        }
        
        c_bar=c_bar/Nbin_max;
        C_bar[t_phys]=c_bar;
        fprintf(file_C_ensemble, "%lf\n", C_bar[t_phys]);
    }

    fclose(file_C_binned);
    fclose(file_C_ensemble);

    /*just a prototype, it will be a function*/
    FILE *file_Energy;
    file_Energy = fopen("../../data_analysis/first_energy_gap.txt", "wt");
    if( file_Energy ==NULL ) {
        perror("Error in opening file");
        exit(1);
    }

    /*funzione che calcola il gap tra gs e primo stato eccitato*/
    double energy_gap [N];

    for (t_phys = 0; t_phys < N; t_phys++)
    {
        double X= (C_bar[(t_phys+N-1)%N]+C_bar[(t_phys+1)%N])/(2*C_bar[t_phys]); 
        energy_gap[t_phys] = log(X + sqrt(X*X-1) );
        fprintf(file_Energy, "%lf\n", energy_gap[t_phys]);
    }
    fclose(file_Energy);

    /*Jackknife clustering*/
    double **C_jack = (double**) malloc(N*sizeof(double*)); 
    for(i=0; i<N; i++)
    {
	    C_jack[i] = (double*) malloc(Nbin_max*sizeof(double));
    }

    double C_mean [N];
    double Var [N];

    FILE *jack_var_E;
    jack_var_E = fopen("../../data_analysis/Jack_variance.txt", "wt");
    if( jack_var_E ==NULL ) {
        perror("Error in opening file");
        exit(1);
    }


    for (t_phys=0; t_phys<N; t_phys++)
    {
        int jack_index;
        double c_mean=0;
        for (t_markov = 0; t_markov < Nbin_max; t_markov++)
        {
                c_mean = c_mean + C[t_phys][t_markov];
        }
        c_mean=c_mean/Nbin_max;

        for (jack_index = 0; jack_index<Nbin_max; jack_index++)
        {
            C_jack[t_phys][jack_index] = c_mean - (1./(double)(Nbin_max-1))*(C[t_phys][jack_index]-c_mean);
        }
        
    }
    
    /*calcolo la matrice energia per ogni cluster jackknife*/
    double **E_jack = (double**) malloc(N*sizeof(double*)); 
    for(i=0; i<N; i++)
    {
	    E_jack[i] = (double*) malloc(Nbin_max*sizeof(double));
    }
    
    

    /*funzione che calcola l'energia di jackknife*/
    for ( t_markov = 0; t_markov < Nbin_max; t_markov++) {

        for ( t_phys = 0; t_phys < N; t_phys++) {

           double X = (C_jack[(t_phys+1)%N][t_markov]+C_jack[(t_phys+N-1)%N][t_markov])/(2*C_jack[t_phys][t_markov]);
           E_jack[t_phys][t_markov] = log(X + sqrt(X*X-1) );
           

        }
        
    }
    
    /*calcolo la varianza dell'energia per ogni tempo fisico */

    for (t_phys=0; t_phys<N;t_phys++){
        double var=0;
        
        for (t_markov = 0; t_markov < Nbin_max; t_markov++) {
            /*calcolo la varianza del cluster jackknife al tempo fisico fissato*/
            var = var + (E_jack[t_phys][t_markov]-energy_gap[t_phys])*(E_jack[t_phys][t_markov]-energy_gap[t_phys]);
        }

        var = var * ((double)Nbin_max-1)/(double)Nbin_max;
        Var[t_phys]=var;
        fprintf(jack_var_E, "%lf\n", var);
    }

    fclose(jack_var_E);

    double Energy=0;
    double weights=0;
    /*calcolo la miglior stima per l'estimatore energia */
    for (t_phys=1; t_phys<9;t_phys++)
    {
        Energy = Energy + ((1/Var[t_phys])*energy_gap[t_phys]);
        weights = weights + 1/Var[t_phys];
    }
    Energy = Energy/weights;
    printf("%lf\n", Energy);

    double *Energy_cluster = (double*) malloc(Nbin_max*sizeof(double*)); 
    
    
    for (t_markov=0; t_markov<Nbin_max; t_markov++)
    {
        double e_cluster =0;
        weights = 0;
        for (t_phys=1; t_phys<9; t_phys++)
        {
            e_cluster = e_cluster + E_jack[t_phys][t_markov] * (1./Var[t_phys]);
            weights = weights + 1./Var[t_phys];
        }

        Energy_cluster[t_markov]=e_cluster/weights;
    }
    
    double variance = 0;

    for (t_markov=1; t_markov<9; t_markov++)
    {
        variance = variance + ((Energy_cluster[t_markov]-Energy)*(Energy_cluster[t_markov]-Energy));
        
    }
    variance = variance *((double)Nbin_max-1.)/(double)Nbin_max;
    printf("%lf\n", sqrt(variance));

    /*funzione che calcola la media per ogni tempo fisico pesando per la relativa varianza il valore di C*/

    return 0;
}
