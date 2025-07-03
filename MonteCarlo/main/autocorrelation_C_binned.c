
#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/global.h"
#include "../include/sweep.h"
#include "../include/random.h"
#include "../include/dazione.h"
#include "../include/azione.h"


int main(){

  int Nbin_max=100000;
  int t_markov; 
  int t_phys;
  double value;

  int cnfg;

  double **BINNED_CORRELATION = (double**) calloc(Nbin_max,sizeof(double*));
  double *C_MEAN = (double*) calloc(N,sizeof(double*));
 
  FILE *correlation;
  correlation = fopen("../../data_analysis/data/harmonic/64/corr_bin_64.txt", "r");
  if (correlation == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }

  FILE *autocorrelation_binned;
  autocorrelation_binned = fopen("../../data_analysis/binned_Gamma.txt", "wt");
  if (autocorrelation_binned == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
    
    /*
    *****************************************************************************
    *
    * Leggiamo i valori del correlatore dal file corr_bin.txt, generati
    * nel main program CORR_FUN_C.
    * Ogni Dbin iterazioni della funzione fscanf viene aggiornato l'indice di riga 
    * Nbin (tempo markoviano) della matrice di correlazione binnata.
    * Leggiamo la varianza del correlatore e la correlazione media dai file 
    * necessarie per il Clustering Jackknife, e calcoliamo la media del gap dE
    * a partire dal valor medio del correlatore per ogni tempo fisico.
    *
    *****************************************************************************
    */

  for (t_markov=0; t_markov<Nbin_max; t_markov++)
  {
    BINNED_CORRELATION[t_markov] = (double*) calloc(N,sizeof(double));
    for (t_phys = 0; t_phys < N; t_phys++)
    {
      fscanf(correlation, "%lf\n", &value);
      BINNED_CORRELATION[t_markov][t_phys] = value;

    }
  }

  /*debug: provo a implementare la media a partire dal correlatore binnato
    non cambia nulla! */
  for (t_phys = 0; t_phys < N; t_phys++)
  {
    double c = 0;
    for (t_markov=0; t_markov<Nbin_max; t_markov++)
    {
      c = c + BINNED_CORRELATION[t_markov][t_phys];
    }
    C_MEAN[t_phys]= c/(double)Nbin_max;
  }


  /*AUTOCORRELATION for the binned correlator*/
  for (t_phys =1; t_phys < 2; t_phys++) 
  {
    /*nota che gamma (T_p = 1 e T_p=63 devono essere uguali: N-i= i )*/
    double N_s_tilde=7000;             /*N_s_tilde deve essere più piccolo di Nbin_max, in particolare */
                              
    for (t_markov = 0; t_markov < N_s_tilde; t_markov++) 
    {

      double g = 0;
      double g_0 =0;

      for (cnfg = 0; cnfg < Nbin_max-N_s_tilde; cnfg++)
        {
          g  = g + (BINNED_CORRELATION[cnfg][t_phys]*BINNED_CORRELATION[cnfg+t_markov][t_phys]);
          g_0 = g_0 +(BINNED_CORRELATION[cnfg][t_phys]*BINNED_CORRELATION[cnfg][t_phys]);
        }
        
      g = g/(Nbin_max-N_s_tilde)-(C_MEAN[t_phys]*C_MEAN[t_phys]);
      g_0 = g_0/(Nbin_max-N_s_tilde)-(C_MEAN[t_phys]*C_MEAN[t_phys]);

      fprintf(autocorrelation_binned, "%d ", t_markov);
      fprintf(autocorrelation_binned, "%lf\n", g/g_0);
    }

    fclose(correlation);
    fclose(autocorrelation_binned);


  }
}
