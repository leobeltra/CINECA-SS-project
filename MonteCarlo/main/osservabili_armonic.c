/*
############################  @ OBSERVABLES_HARMONIC_C @   ################################
#
# Federico De Matteis, Lab. di Fisica Computazionale 2022/2023
# Calcolo delle osservabili fisiche e del loro errore (dE_01 e |<0|x|1>|^2)
# per il sistema armonico.
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

void free_fun(double **MATRIX, int raws)
{
  int raw;
  for (raw = 0; raw < raws; raw++)
  {
    free(MATRIX[raw]);
  }
  free(MATRIX);
}


int main(){

  int Dbin=50;                                                                /*larghezza di un bin nella procedura di binning (in corr_fun.c)*/
  int Nbin_max=20000;                                                         /*numero di configurazioni indipendenti generate dal Metropolis (in corr_fun.c)*/
  int t_markov;                                                                /*tempo Markoviano*/
  int t_phys;                                                                  /*distanza delle sorgenti nel tempo fisico*/

  double value;                                                                /*variabili temporanee*/
  double value1, value2;
  double value_3;
  int Nbin=0;
  int cnfg;                                                                  
  double PHYSICAL_ENERGY=0;       
  double weights_E=0;
  double weights_M=0;
  double PHYSICAL_VARIANCE_E = 0;
  double PHYSICAL_VARIANCE_M = 0;
  double PHYSICAL_M=0;
  
  int sanitizer_E=0;                                                            /*variabili per sanificazione delle misure*/
  int sanitizer_M=0;


  double **BINNED_CORRELATION = (double**) calloc(Nbin_max,sizeof(double*));    /*matrice del correlatore binnato*/
  double **BINNED_CORRELATION_2 = (double**) calloc(Nbin_max,sizeof(double*));
  double **C_JACK = (double**) calloc(Nbin_max,sizeof(double*));                /*cluster jackknife del correlatore*/
  double **E_JACK = (double**) calloc(Nbin_max,sizeof(double*));                /*cluster jackknife per il gap di energia */
  double **MATRIX_EL_JACK = (double**) calloc(Nbin_max,sizeof(double*));        /*cluster jackknife per l'elemento di matrice*/

  double *C_MEAN = (double*) calloc(N,sizeof(double*));                         /*correlatore medio*/
  double *VAR_C = (double*) calloc(N,sizeof(double*));                          /*varianza del correlatore medio*/
  double *DELTA_E = (double*) calloc(N,sizeof(double*));                        /*media pesata del gap di energia (migliore stima)*/
  double *VAR_dE = (double*) calloc(N,sizeof(double*));                         /*varianza del gap di energia*/
  double *VAR_M = (double*) calloc(N,sizeof(double*));                          /*media pesata dell'elemento di matrice (migliore stima)*/
  double *MATRIX_EL = (double*) calloc(N,sizeof(double*));                      /*varianza dell'elemento di matrice*/

  double *Energy_cluster = (double*) malloc(Nbin_max*sizeof(double*));          /*cluster jackknife per il gap di energia*/
  double *Matrix_cluster = (double*) malloc(Nbin_max*sizeof(double*));          /*cluster jackknife per l'elemento di matrice*/
  

  FILE *correlation;
  correlation = fopen("../../data_analysis/data/harmonic/64/corr_bin_64.txt", "r");
  if (correlation == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
  /*implementare la funzione di correlazione a due punti <x_l^2 x_k^2>*/
  /*FILE *correlation_2;
  correlation_2 = fopen("../../data_analysis/data/harmonic/64/corr_bin_64.txt", "r");
  if (correlation == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }*/
  FILE *corr_mean;
  corr_mean = fopen("../../data_analysis/data/harmonic/64/corr_mean_64.txt", "r");
  if (corr_mean == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
  FILE *variance;
  variance = fopen("../../data_analysis/data/harmonic/64/variance_64.txt", "r");
  if (variance == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
  
  FILE *file_energy;
  file_energy = fopen("../../data_analysis/energy_measure_64.txt", "wt");
  if (file_energy == NULL)
  {
    printf("Problem in opening file in wt mode.\n");
    exit(1);
  }
  FILE *file_matrix;
  file_matrix = fopen("../../data_analysis/matrix_measure_64.txt", "wt");
  if (file_matrix == NULL)
  {
    printf("Problem in opening file in wt mode.\n");
    exit(1);
  }
  FILE *jack_var_dE;
  jack_var_dE = fopen("../../data_analysis/jack_var_dE_64.txt", "wt");
  if (file_energy == NULL)
  {
    printf("Problem in opening file in wt mode.\n");
    exit(1);
  }
  FILE *jack_var_M;
  jack_var_M = fopen("../../data_analysis/jack_var_M_64.txt", "wt");
  if (jack_var_M == NULL)
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
    
  FILE *sigma_C;
  sigma_C = fopen("../../data_analysis/sigma_C_medio_64.txt", "wt");
  if (sigma_C == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
  
  FILE *sigma_relative;
  sigma_relative = fopen("../../data_analysis/sigma_C_relative_64.txt", "wt");
  if (sigma_relative == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
  
  FILE *c_medio;
  c_medio = fopen("../../data_analysis/c_medio_64.txt", "wt");
  if (c_medio == NULL)
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
    /*BINNED_CORRELATION_2[t_markov] = (double*) calloc(N,sizeof(double));  */  
    for (t_phys = 0; t_phys < N; t_phys++)
    {
      fscanf(correlation, "%lf\n", &value);
      BINNED_CORRELATION[t_markov][t_phys] = value;
      
      /*fscanf(correlation_2, "%lf\n", &value_3);      
      BINNED_CORRELATION_2[t_markov][t_phys] = value_3;*/
    }
  }

  /*debug: provo a implementare la media a partire dal correlatore binnato
    non cambia nulla! :() */
  for (t_phys = 0; t_phys < N; t_phys++)
  {

    double c = 0;
    for (t_markov=0; t_markov<Nbin_max; t_markov++)
    {
      c = c + BINNED_CORRELATION[t_markov][t_phys];

    }
    C_MEAN[t_phys]= c/(double)Nbin_max;
    /*printf("%lf\n", C_MEAN[t_phys]);*/
  }
  
  /*calcolo dell'errore di <C(t)> e dell'errore relativo*/
  for (t_phys = 0; t_phys < N; t_phys++)
  {
    double mean_square = 0;

    for (t_markov=0; t_markov<Nbin_max; t_markov++)
    {
      mean_square = mean_square+(BINNED_CORRELATION[t_markov][t_phys]-C_MEAN[t_phys])*(BINNED_CORRELATION[t_markov][t_phys]-C_MEAN[t_phys]);
    }

    double var_c_mean= mean_square/(double)((Nbin_max-1));

    fprintf(c_medio, "%d %lf %lf\n", t_phys, C_MEAN[t_phys], sqrt(var_c_mean/(double)Nbin_max)); 
    fprintf(sigma_C, "%d %lf\n", t_phys, sqrt(var_c_mean/(double)Nbin_max));                         
    fprintf(sigma_relative, "%d %lf\n", t_phys, sqrt(var_c_mean/(double)Nbin_max)/C_MEAN[t_phys]);    
  }

  /*
  for (t_phys = 0; t_phys < N; t_phys++)
  {
    fscanf(variance, "%lf\n", &value1);
    VAR_C[t_phys]=(value1*value1)/(double)Nbin_max; 

    fprintf(c_medio, "%d %lf %lf\n", t_phys, C_MEAN[t_phys], sqrt(VAR_C[t_phys])); 
    fprintf(sigma_C, "%d %lf\n", t_phys, sqrt(VAR_C[t_phys]));                          
    fprintf(sigma_relative, "%d %lf\n", t_phys, sqrt(VAR_C[t_phys])/C_MEAN[t_phys]);   
  }*/

/*
*****************************************************************************
*
* Creazione dei cluster Jackknife per la matrice del correlatore calcolo
* della media del gap energetico e del quadrato dell'elemento di matrice.
*
*****************************************************************************
*/

  for (t_markov=0; t_markov<Nbin_max; t_markov++)
  {
    C_JACK[t_markov] = (double*) malloc(N*sizeof(double));
  }
 
  for (t_phys=0; t_phys<N; t_phys++)
  {
    int jack_index;
    double c_mean=0;
    
    double X = (C_MEAN[(t_phys+N-1)%N]+C_MEAN[(t_phys+1)%N])/(2*C_MEAN[t_phys]); /*si implementa Eq() in relazione*/
    DELTA_E[t_phys] = acosh(X);/*log(X + sqrt(X*X-1) );*/                        /*per il calcolo dell'energia*/
    /*fprintf(file_energy, "%d ", t_phys);
    fprintf(file_energy, "%lf\n", DELTA_E[t_phys]/SPACING);*/

    double Z = cosh( ((double)N/2.-t_phys) * DELTA_E[t_phys] );                  /*si implementa Eq() in relazione*/
    double Y = C_MEAN[t_phys] * exp(((double)N/2)*DELTA_E[t_phys]);              /*per il calcolo dell'elemento di matrice al quadrato*/
    
    MATRIX_EL[t_phys] = Y/(2*Z);
    /*fprintf(file_matrix, "%d ", t_phys);
    fprintf(file_matrix, "%lf\n", MATRIX_EL[t_phys]);*/                          /*qui comincia l'analisi degli errori*/
    for (jack_index = 0; jack_index<Nbin_max; jack_index++)                      /*creo le variabili di Jackknife primarie C^j(t) dove t è il tempo fisico*/
    {
      C_JACK[jack_index][t_phys] = C_MEAN[t_phys] - (1./(double)(Nbin_max-1))*(BINNED_CORRELATION[jack_index][t_phys]-C_MEAN[t_phys]);
    }
  }
  
  for ( t_markov = 0; t_markov < Nbin_max; t_markov++) 
  {
    E_JACK[t_markov] = (double*) malloc(N*sizeof(double));
    MATRIX_EL_JACK[t_markov] = (double*) malloc(N*sizeof(double));

    for ( t_phys = 0; t_phys < N; t_phys++) {                                  /* si implementa  Equazione () in relazione*/
                                                                    
      double X = (C_JACK[t_markov][(t_phys+1)%N]+C_JACK[t_markov][(t_phys+N-1)%N])/(2*C_JACK[t_markov][t_phys]); 
      E_JACK[t_markov][t_phys] = acosh(X);/*log(X + sqrt(X*X-1) );*/

      double Z = cosh(((double)N/2-t_phys)*E_JACK[t_markov][t_phys]);         /*si implementa   Equazione () in relazione*/
      double Y = C_JACK[t_markov][t_phys] * exp(((double)N/2)*E_JACK[t_markov][t_phys]);
      
      MATRIX_EL_JACK[t_markov][t_phys] = Y/(2*Z);

    }
  }
  
  free_fun(BINNED_CORRELATION, Nbin_max);

  /*
  *****************************************************************************
  *
  * La memoria occupata dal correlatore binnato è stata liberata dopo aver
  * creato i cluster jackknife per il correlatore binnato e per  <E_1|dE|E_0> 
  * e <E_1|x|E_0>; 
  * Passiamo a calcolare la varianza di <E_1|dE|E_0> e <E_1|x|E_0> a partire 
  * dai cluster jackknife.
  * 
  *****************************************************************************
  */

  for (t_phys=0; t_phys<N;t_phys++)
  {
    double var_E=0;
    double var_M=0;
    double dE= DELTA_E[t_phys]; /*DELTA_E è il gap di E. calcolato dalla correlazione binnata*/
    double M_EL = MATRIX_EL[t_phys];

    for (t_markov = 0; t_markov < Nbin_max; t_markov++) 
    {
      var_E = var_E + (E_JACK[t_markov][t_phys]-dE)*(E_JACK[t_markov][t_phys]-dE);                      /*si implementa Eq () in relazione*/
      var_M = var_M + (MATRIX_EL_JACK[t_markov][t_phys]-M_EL)*(MATRIX_EL_JACK[t_markov][t_phys]-M_EL);  /*si implementa Eq () in relazione*/
    }
    var_E = var_E * ((double)Nbin_max-1.)/(double)Nbin_max;
    var_M = var_M * ((double)Nbin_max-1.)/(double)Nbin_max;
    VAR_dE[t_phys]=var_E;
    VAR_M[t_phys]= var_M; /*questi 2 vettori non sono necessari, salvo gia su file */

    fprintf(jack_var_dE, "%lf\n", sqrt(VAR_dE[t_phys])/SPACING);
    fprintf(jack_var_M, "%lf\n", sqrt(VAR_M[t_phys]));
  }

  free_fun(C_JACK, Nbin_max);

  /*
  *****************************************************************************
  *
  * Per ogni tempo fisico abbiamo quindi una misura di energia con la 
  * relativa varianza, trovata tramite la procedura di Jackknife (abbiamo un vettore di 
  * misure dE di dimensione N).
  * Ora siamo in grado di eseguire una media del gap di energia, pesando ogni 
  * contributo ad un determinato tempo fisico per l'inverso della relativa varianza
  * 
  * Lo stesso procedimento viene usato per l'elemento di matrice al quadrato.
  * 
  *****************************************************************************
  */
  

  for (t_phys=0; t_phys<N/2;t_phys++){
    
    if (!isnan(DELTA_E[t_phys]) && !isnan(VAR_dE[t_phys]) 
    && DELTA_E[t_phys] < 10 && DELTA_E[t_phys] > 0 && VAR_dE[t_phys]>0)
    {
      if(sqrt(VAR_dE[t_phys])/DELTA_E[t_phys] < 0.3) /*scarto le misure non accettabili*/
      {
        PHYSICAL_ENERGY = PHYSICAL_ENERGY + ((1/VAR_dE[t_phys])*DELTA_E[t_phys]);
        weights_E = weights_E + 1./VAR_dE[t_phys];
        sanitizer_E = sanitizer_E +1;

        printf("sigma_dE  | dE\n");
        printf("%lf  |  ", sqrt(VAR_dE[t_phys])/SPACING);
        printf("%lf\n", DELTA_E[t_phys]/SPACING);
      }
    }
    
    if (!isnan(MATRIX_EL[t_phys]) && !isnan(VAR_M[t_phys]) 
    && MATRIX_EL[t_phys] < 10 && MATRIX_EL[t_phys] > 0 && VAR_M[t_phys]>0)
    {
      if(sqrt(VAR_M[t_phys])/MATRIX_EL[t_phys] < 0.3) /*scarto le misure non accettabili*/
      {
        PHYSICAL_M = PHYSICAL_M + ((1/VAR_M[t_phys])*MATRIX_EL[t_phys]);
        weights_M = weights_M + 1./VAR_M[t_phys];
        sanitizer_M=sanitizer_M + 1;

        printf("var_M^2  |  M^2\n");
        printf("%lf  |  ", VAR_M[t_phys]);
        printf("%lf\n", MATRIX_EL[t_phys]);
        
      }
    }
    
  }
  
  if (sanitizer_E == 0)
  {
    printf("Nessuna misura di E_1-E_0 utilizzabile \n");
    printf("Si consiglia di aumentare il sampling, modificando:\n");
    printf("Dbin, Nbin_max, M_sweep(global)\n");
    printf("nota: deve valere Dbin*Nbin_max=M_sweep\n");
  }

  if(sanitizer_M==0){
    printf("Nessuna misura di |M_{01}|^2 utilizzabile.\n");
    printf("Si consiglia di aumentare il sampling, modificando:\n");
    printf("Dbin, Nbin_max, M_sweep(global)\n");
    printf("nota: deve valere Dbin*Nbin_max=M_sweep\n");
  }

  if(sanitizer_E==0 && sanitizer_M==0){
    return 0;
  }

  /* calcolo delle varianze delle osservabili terziarie tramite jackknife*/
  PHYSICAL_ENERGY = PHYSICAL_ENERGY/weights_E;
  PHYSICAL_M = PHYSICAL_M/weights_M;

  for (t_markov=0; t_markov<Nbin_max; t_markov++)
  {
    double e_cluster =0;
    double m_cluster =0;

    weights_E = 0;
    weights_M = 0;
    
    for (t_phys=0; t_phys<N/2; t_phys++)
    {
      if(!isnan(DELTA_E[t_phys]) && !isnan(VAR_dE[t_phys]) 
      && sqrt(VAR_dE[t_phys])/DELTA_E[t_phys] < 0.3 
      && DELTA_E[t_phys] < 10 && DELTA_E[t_phys] > 0 && VAR_dE[t_phys]>0)
      {
        e_cluster = e_cluster + E_JACK[t_markov][t_phys] * (1./VAR_dE[t_phys]);
        weights_E = weights_E + 1./VAR_dE[t_phys];
      }
      

      if(!isnan(MATRIX_EL[t_phys]) && !isnan(VAR_M[t_phys]) 
      && sqrt(VAR_M[t_phys])/MATRIX_EL[t_phys] < 0.3 
      && MATRIX_EL[t_phys]<10  && MATRIX_EL[t_phys]> 0 && VAR_M[t_phys]>0)
      {
        m_cluster= m_cluster+MATRIX_EL_JACK[t_markov][t_phys]*(1./VAR_M[t_phys]);
        weights_M = weights_M + 1./VAR_M[t_phys];
        /*printf("varianze usate nel secondo ciclo \n");
        printf("%lf\n", VAR_M[t_phys]);
        printf("elementi di matrice usati nel secondo ciclo\n");
        printf("%lf\n", MATRIX_EL[t_phys]);*/
      } 
    
    }

    /*++++++++++++++++*/

    Energy_cluster[t_markov]=e_cluster/weights_E;
    Matrix_cluster[t_markov]=m_cluster/weights_M;
  }

  /*for (t_markov=0; t_markov<Nbin_max; t_markov++)
  {
    printf("%lf\n", Energy_cluster[t_markov]);
  }*/
  /*printf("physical E");
  printf("%lf\n", PHYSICAL_ENERGY);*/
  

  free_fun(E_JACK, Nbin_max);
  free_fun(MATRIX_EL_JACK, Nbin_max);

  for (t_markov=0; t_markov<Nbin_max; t_markov++)
  {
    PHYSICAL_VARIANCE_E = PHYSICAL_VARIANCE_E + ((Energy_cluster[t_markov]-PHYSICAL_ENERGY)*(Energy_cluster[t_markov]-PHYSICAL_ENERGY));
    PHYSICAL_VARIANCE_M = PHYSICAL_VARIANCE_M +(Matrix_cluster[t_markov]-PHYSICAL_M)*(Matrix_cluster[t_markov]-PHYSICAL_M);
  }

  PHYSICAL_VARIANCE_E = PHYSICAL_VARIANCE_E*(((double)Nbin_max-1.)/(double)Nbin_max);
  PHYSICAL_VARIANCE_M = PHYSICAL_VARIANCE_M*(((double)Nbin_max-1.)/(double)Nbin_max);


  /*
  *****************************************************************************
  * 
  * Libero la memoria utilizzata residua.
  * 
  *****************************************************************************
  */
  for (t_phys=0; t_phys<N; t_phys++)

  {
    fprintf(file_energy, "%d ", t_phys);
    fprintf(file_energy, "%lf ", DELTA_E[t_phys]/SPACING);
    fprintf(file_energy, "%lf\n", sqrt(VAR_dE[t_phys])/SPACING);
    fprintf(file_matrix, "%d ", t_phys);
    fprintf(file_matrix, "%lf ", MATRIX_EL[t_phys]);
    fprintf(file_matrix, "%lf\n", sqrt(VAR_M[t_phys]));
  }
  

  free(Energy_cluster);
  free(C_MEAN);
  free(VAR_C);
  free(VAR_dE);
  free(DELTA_E);
  fclose(sigma_C);
  fclose(c_medio);
  fclose(sigma_relative);
  fclose(correlation);
  fclose(corr_mean);
  fclose(variance);
  fclose(file_energy);
  fclose(jack_var_dE);
  fclose(jack_var_M);
  fclose(file_matrix);

  /*fclose(correlation_2);*/
  /*fclose(autocorrelation_binned);*/

  /*calcolo dei valori teorici nel discreto di dE_{01} e M_{01}^2 
    e confroto con i valori simulati tramite la variabile t-student*/

  double omega_bar = sqrt((1+SPACING*SPACING/4));                                                             /*risultato teorico in Eq() in relazione*/
  double omega_tilde = (1./SPACING)*log(SPACING*sqrt(1+(SPACING*SPACING)/4.)+1.+(SPACING*SPACING/2.));        /*risultato teorico in Eq() in relazione*/
  double M_th_lattice_squared = (1./sqrt(4+SPACING*SPACING));                                                 /*risultato teorico in Eq() in relazione*/

  printf("valore teorico di M su reticolo");
  printf("%lf\n", M_th_lattice_squared);
  printf("valore teorico di dE su reticolo");
  printf("%lf\n", omega_tilde);

  double t_student_E = ((PHYSICAL_ENERGY/SPACING) - omega_tilde);
  t_student_E=fabs(t_student_E)/(sqrt(PHYSICAL_VARIANCE_E)/SPACING);         /*calcolo delle t-student*/
  printf("%lf\n", t_student_E);

  double t_student_M = fabs(PHYSICAL_M - M_th_lattice_squared)/sqrt(PHYSICAL_VARIANCE_M);
  
  printf("*******************************************************\n");

  printf("La media del gap di energia è ");
  printf("%lf\n", PHYSICAL_ENERGY/SPACING);

  printf("La media dell'elemento di matrice è ");
  printf("%lf\n", PHYSICAL_M);

  printf("L'errore del gap di energia è ");
  printf("%lf\n", sqrt(PHYSICAL_VARIANCE_E)/SPACING);

  printf("L'errore dell'elemento di matrice è ");
  printf("%lf\n", sqrt(PHYSICAL_VARIANCE_M));

  printf("*******************************************************\n");
  printf("t-student dE ");
  printf("%lf\n", t_student_E);

  printf("t-student M ");
  printf("%lf\n", t_student_M);



  return 0;
}
