/*
############################  @ OSSERVABILI_C @   ################################
#
# Federico De Matteis, Lab. di Fisica Computazionale 2022/2023
# Calcolo delle osservabili fisiche e del loro errore (dE_01 e |<0|x|1>|^2) per 
# il sistema anarmonico
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

  int Dbin=150;
  int Nbin_max=500000;
  int t_markov; 
  int t_phys;
  double value;
  double value1, value2;
  int Nbin=0;
  double PHYSICAL_ENERGY=0;
  double weights_E=0;
  double weights_M=0;
  double PHYSICAL_VARIANCE_E = 0;
  double PHYSICAL_VARIANCE_M = 0;

  double **BINNED_CORRELATION = (double**) calloc(Nbin_max,sizeof(double*));
  double **C_JACK = (double**) calloc(Nbin_max,sizeof(double*));
  double **E_JACK = (double**) calloc(Nbin_max,sizeof(double*));
  double **MATRIX_EL_JACK = (double**) calloc(Nbin_max,sizeof(double*));

  double *C_MEAN = (double*) calloc(N,sizeof(double*));
  /*in realtà la variaza del correlatore non mi serve*/
  double *VAR_C = (double*) calloc(N,sizeof(double*));
  /*varianza di dE*/
  double *DELTA_E = (double*) calloc(N,sizeof(double*));
  double *VAR_dE = (double*) calloc(N,sizeof(double*));
  double *VAR_M = (double*) calloc(N,sizeof(double*));
  double *MATRIX_EL = (double*) calloc(N,sizeof(double*));  

  double *Energy_cluster = (double*) malloc(Nbin_max*sizeof(double*));
  double *Matrix_cluster = (double*) malloc(Nbin_max*sizeof(double*));
  

  FILE *correlation;
  correlation = fopen("../../data_analysis/data/harmonic/128/corr_bin_128_final.txt", "r");
  if (correlation == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
  FILE *corr_mean;
  corr_mean = fopen("../../data_analysis/data/harmonic/128/corr_mean_128_final.txt", "r");
  if (corr_mean == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
  FILE *variance;
  variance = fopen("../../data_analysis/data/harmonic/128/variance_128_final.txt", "r");
  if (variance == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
  FILE *file_energy;
  file_energy = fopen("../../data_analysis/energy_measure.txt", "wt");
  if (file_energy == NULL)
  {
    printf("Problem in opening file in wt mode.\n");
    exit(1);
  }
  FILE *file_matrix;
  file_matrix = fopen("../../data_analysis/matrix_measure.txt", "wt");
  if (file_matrix == NULL)
  {
    printf("Problem in opening file in wt mode.\n");
    exit(1);
  }
  FILE *jack_var_dE;
  jack_var_dE = fopen("../../data_analysis/jack_var_dE.txt", "wt");
  if (file_energy == NULL)
  {
    printf("Problem in opening file in wt mode.\n");
    exit(1);
  }
  FILE *jack_var_M;
  jack_var_M = fopen("../../data_analysis/jack_var_M.txt", "wt");
  if (jack_var_M == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
  
  FILE *sigma_c_mean;
  sigma_c_mean = fopen("../../data_analysis/sigma_c_mean.txt", "wt");
  if (sigma_c_mean == NULL)
  {
    printf("Cannot open file for reading.\n");
    exit(1);
  }
/*
*****************************************************************************
*
* Leggiamo i valori del correlatore dal file corr_bin.txt, generati
* nel main program CORR_FUN_C.
* Ogni N iterazioni della funzione fscanf viene aggiornato l'indice di riga 
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

  for (t_phys = 0; t_phys < N; t_phys++)
  {
    fscanf(variance, "%lf\n", &value1);
    VAR_C[t_phys]=value1/sqrt((double)Nbin_max);
    /*fscanf(corr_mean, "%lf\n", &value2);
    C_MEAN[t_phys]=value2;*/
    fprintf(sigma_c_mean,"%lf\n", VAR_C[t_phys]);
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
  }
  

/*
*****************************************************************************
*
* Creazione dei cluster Jackknife per la matrice del correlatore calcolo
* della media delgap energetico e del quadrato dell'elemento di matrice.
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
    
    double X= (C_MEAN[(t_phys+N-1)%N]+C_MEAN[(t_phys+1)%N])/(2*C_MEAN[t_phys]); 
    DELTA_E[t_phys] = acosh(X);/*log(X + sqrt(X*X-1) );*/
    fprintf(file_energy, "%lf\n", DELTA_E[t_phys]/SPACING);

    double Z = cosh( ((double)N/2.-t_phys) * DELTA_E[t_phys] );
    double Y = C_MEAN[t_phys] * exp(((double)N/2)*DELTA_E[t_phys]);
    
    MATRIX_EL[t_phys] = Y/(2*Z);
    fprintf(file_matrix, "%lf\n", MATRIX_EL[t_phys]);

    for (jack_index = 0; jack_index<Nbin_max; jack_index++)
    {
      C_JACK[jack_index][t_phys] = C_MEAN[t_phys] - (1./(double)(Nbin_max-1))*(BINNED_CORRELATION[jack_index][t_phys]-C_MEAN[t_phys]);
    }
  }
  
  for ( t_markov = 0; t_markov < Nbin_max; t_markov++) 
  {
    E_JACK[t_markov] = (double*) malloc(N*sizeof(double));
    MATRIX_EL_JACK[t_markov] = (double*) malloc(N*sizeof(double));

    for ( t_phys = 0; t_phys < N; t_phys++) {

      double X = (C_JACK[t_markov][(t_phys+1)%N]+C_JACK[t_markov][(t_phys+N-1)%N])/(2*C_JACK[t_markov][t_phys]);
      E_JACK[t_markov][t_phys] = acosh(X);/*log(X + sqrt(X*X-1) );*/

      double Z = cosh(((double)N/2-t_phys)*E_JACK[t_markov][t_phys]);
      double Y = C_JACK[t_markov][t_phys] * exp(((double)N/2)*E_JACK[t_markov][t_phys]);
      
      MATRIX_EL_JACK[t_markov][t_phys] = Y/(2*Z);

    }
  }
  
  free_fun(BINNED_CORRELATION, Nbin_max);

  /*
  *****************************************************************************
  *
  * Abbiamo liberato la memoria occupata dal correlatore binnato dopo aver
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
      var_E = var_E + (E_JACK[t_markov][t_phys]-dE)*(E_JACK[t_markov][t_phys]-dE);
      var_M = var_M + (MATRIX_EL_JACK[t_markov][t_phys]-M_EL)*(MATRIX_EL_JACK[t_markov][t_phys]-M_EL);
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
  * Per ogni tempo fisico abbiamo quindi  una misura di energia con la 
  * relativa varianza, trovata tramite la procedura di Jackknife.
  * Ora siamo in grado di eseguire una media del gap di energia, pesando ogni 
  * contributo ad un determinato tempo fisico per l'inverso della relativa varianza; 
  * Assumeremo il valore trovato come il coefficiente della la retta interpolante
  * che rappresenta la nostra migliore stima di un estimatore maximum likelihood 
  * del gap di energia tra lo stato di vuoto e il primo stato eccitato.
  * Lo stesso procedimento viene usato per l'elemento di matrice al quadrato.
  * 
  *****************************************************************************
  */
  
  
  double PHYSICAL_M=0;
  /*questo for viene sostituito dal successivo*/
  /*for (t_phys=1; t_phys<4;t_phys++)
  {
    PHYSICAL_ENERGY = PHYSICAL_ENERGY + ((1/VAR_dE[t_phys])*DELTA_E[t_phys]);
    weights_E = weights_E + 1./VAR_dE[t_phys];

    PHYSICAL_M = PHYSICAL_M + ((1/VAR_M[t_phys])*MATRIX_EL[t_phys]);
    weights_M = weights_M+1./VAR_M[t_phys];
  }
  /***************************************/
  int sanitizer_E=0;
  int sanitizer_M=0;

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

  /*++++++++++++++++ calcolo delle varianze delle osservabili a partire dal jackknife*/
  PHYSICAL_ENERGY = PHYSICAL_ENERGY/weights_E;
  PHYSICAL_M = PHYSICAL_M/weights_M;

  for (t_markov=0; t_markov<Nbin_max; t_markov++)
  {
    double e_cluster =0;
    double m_cluster =0;

    weights_E = 0;
    weights_M = 0;
    
    /*questo for viene sostituito dal successivo*/
    /*for (t_phys=1; t_phys<4; t_phys++)
    {
      e_cluster = e_cluster + E_JACK[t_markov][t_phys] * (1./VAR_dE[t_phys]);
      weights_E = weights_E + 1./VAR_dE[t_phys];
      m_cluster= m_cluster+MATRIX_EL_JACK[t_markov][t_phys]*(1./VAR_M[t_phys]);
      weights_M = weights_M + 1./VAR_M[t_phys];
    }
    /*++++++++++++++++*/

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
  free(Energy_cluster);
  free(C_MEAN);
  free(VAR_C);
  free(VAR_dE);
  free(DELTA_E);

  fclose(correlation);
  fclose(corr_mean);
  fclose(variance);
  fclose(file_energy);
  fclose(jack_var_dE);
  fclose(jack_var_M);
  fclose(file_matrix);
  fclose(sigma_c_mean);

  /*calcolo i valori teorici di dE e |M_{01}|^2 dalla teoria 
    delle perturbazioni degenere, per l'oscillatore anarmonico;
    Calcolo la t-student per gni misura, rispetto al valore vero*/

    double dE_an_th_lattice_1 = (1./SPACING)*log(SPACING*sqrt(1+(SPACING*SPACING)/4.)+1.+(SPACING*SPACING/2.))
    +3.*lambda/(1.+SPACING*SPACING/4.);

    double dE_an_th_lattice_2 = (1./SPACING)*log(SPACING*sqrt(1+(SPACING*SPACING)/4.)+1.+(SPACING*SPACING/2.))
    +3.*lambda/(1.+SPACING*SPACING/4.)
    -18.*lambda*lambda/((1.+SPACING*SPACING/4)*(1.+SPACING*SPACING/4)* ((1./SPACING)*log(SPACING*sqrt(1+(SPACING*SPACING)/4.)+1.+(SPACING*SPACING/2.))) ) ;
    
    double dE_an_th = 1. + 3*lambda - 18.*lambda*lambda; /*+ (1791/8)*(lambda*lambda*lambda);*/
    
    double t_student_E_1 = ((PHYSICAL_ENERGY/SPACING) - dE_an_th_lattice_1)/(sqrt(PHYSICAL_VARIANCE_E/SPACING));
    double t_student_E_2 = ((PHYSICAL_ENERGY/SPACING) - dE_an_th_lattice_2)/(sqrt(PHYSICAL_VARIANCE_E)/SPACING) ;


    /*calcoli corretti*/
    double M_an_th_lattice = (1./(sqrt(2*sqrt(1+SPACING*SPACING/4))))*(1-(sqrt(2)/(2*(1./SPACING)*log(SPACING*sqrt(1+(SPACING*SPACING)/4.)+1.+(SPACING*SPACING/2.))))*(lambda/(1+SPACING*SPACING/4))*((3+3*sqrt(2))/4) );
    M_an_th_lattice = M_an_th_lattice*M_an_th_lattice;
    double M_an_th = 1/sqrt(2) - lambda/2 * ((3+3*sqrt(2))/4);
    M_an_th = M_an_th*M_an_th;


    double omega_bar = sqrt((1+SPACING*SPACING/4));
    double omega_tilde = (1./SPACING)*log(SPACING*sqrt(1+(SPACING*SPACING)/4.)+1.+(SPACING*SPACING/2.));

    double M_an_lattice_2 = 1./(2*omega_bar) - 3.*lambda/(2*omega_bar*omega_bar*omega_bar*omega_tilde) + 525*lambda*lambda/(32*omega_bar*omega_bar*omega_bar*omega_bar*omega_bar*omega_tilde*omega_tilde);

    /*calcoli a mano*/
    double M_an_th_2 = (1.-lambda*lambda*85/16)*(1./sqrt(2))*(1.-(lambda/(1.))*(sqrt(2)*(3+3*sqrt(2))/4.)/(2.)+(lambda*lambda/(sqrt(2.)))*(75*sqrt(2)/16. -39./64.))
                                +((lambda/(1.))*(5*sqrt(6)/(2.))+(lambda*lambda*sqrt(6))*(195./8.)/2.)*((1./sqrt(2.))*(lambda/(1.))*(((3+3*sqrt(2))/4.)*sqrt(3)/(2.)-(((2*sqrt(3)/4.)*sqrt(4))/(4.))) 
                                + (lambda*lambda/(sqrt(2.)))*((75*sqrt(3)/16)+(18*sqrt(4)/16)) )
                                + (-lambda*(sqrt(30)/2.)/(4.) + lambda*lambda*(33./16.)*sqrt(120))*((-1./sqrt(2.))*(lambda*lambda/(1.))*(2*sqrt(3)*sqrt(5)/(16.)) + (lambda*lambda/sqrt(2.))*(18*sqrt(5)/16 + 17*sqrt(6)/192.) )
                                + ((lambda*lambda/sqrt(2.))*(17*sqrt(7)/192 + sqrt(8)/152)*(23*sqrt(5040)/192))
                                + (1904.940944/512.)*(lambda*lambda*sqrt(9)/(512.*sqrt(2.))); 
    
    M_an_th_2 = M_an_th_2* M_an_th_2;
    
    double t_student_M_1 = (PHYSICAL_M-M_an_th_lattice)/sqrt(PHYSICAL_VARIANCE_M);
    double t_student_M_2 = (PHYSICAL_M-M_an_lattice_2)/sqrt(PHYSICAL_VARIANCE_M);

    printf("*******************************************************\n");
    printf("il valore teorico di dE al primo ordine sul reticolo è ");
    printf("%lf\n", dE_an_th_lattice_1);

    printf("il valore teorico di dE al secondo ordine sul reticolo è ");
    printf("%lf\n", dE_an_th_lattice_2);

    printf("il valore teorico di dE nel continuo è ");
    printf("%lf\n", dE_an_th);

    printf("il valore teorico di |M_01|^2 al primo ordine sul reticolo è ");
    printf("%lf\n", M_an_th_lattice);

    printf("il valore teorico di |M_01|^2 al primo ordine sul reticolo è ");
    printf("%lf\n", M_an_lattice_2);

    printf("il valore teorico di |M_01|^2 nel continuo è ");
    printf("%lf\n", M_an_th_2);


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

  printf("t-student dE al primo ordine è");
  printf("%lf\n", t_student_E_1);

  printf("t-student dE al secondo ordine è");
  printf("%lf\n", t_student_E_2);

  printf("t-student |M_01|^2 al primo ordine è");
  printf("%lf\n", t_student_M_1);

  printf("t-student |M_01|^2 al secondo ordine è");
  printf("%lf\n", t_student_M_2);

  printf("*******************************************************\n");
  printf("T-student  limite al continuo\n");
  double x = 0;
  double f_E = 1.02741 - 0.0435071*x*x + 0.000265915*x*x*x*x;
  double error_extrapolation_E = 0.000335269 + 0.00335129 + 0.00320892;
  double t_student_continuum_E = (f_E - dE_an_th)/error_extrapolation_E;
  printf("t-student dell'estrapolazione al continuo per E è");
  printf("%lf\n", t_student_continuum_E);

  double f_M = 0.486703 - 0.0644067*x*x + 0.0117488*x*x*x*x;
  double error_extrapolation_M = 0.000241261 + 0.00241161 + 0.00230916  ;
  double t_student_continuum_M = (f_M - M_an_th)/error_extrapolation_M;
  printf("t-student dell'estrapolazione al continuo per M è");
  printf("%lf\n", t_student_continuum_M);

  printf("*******************************************************\n");

  return 0;
}
