/*
############################  @ OBSERVABLES_HARMONIC_C @   ################################
#
# Federico De Matteis, Lab. di Fisica Computazionale 2022/2023
# Stato di vuoto dell'oscillatore armonico monodimensionale.
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
#include "../include/dazione_an_creutz.h"

int hash_histo(t_phys)
{
  int Bin;
  Bin = ((int)(100*(xx[t_phys]+2.5)/(5.))-1);
  return Bin;
}

int main()
{
    int HISTO[100];
    int i; 
    for (i = 0; i < 100; i++)
    {
        HISTO[i] = 0;
    }
    
    int t_phys;
    int t_markov;
    double r3[N];
    double r4[N];
    double x_min=-2.5;
    double x_max=2.5;
    
    int Nbin =0;
    int Nbin_max = 3000000;
    int Dbin=40;
    int bin;
    double shift=0;

    FILE *Histo_Q;
    Histo_Q = fopen("../../data_analysis/Histo_Q.txt", "wt");
    if (Histo_Q == NULL)
    {
        printf("Problem in opening file.\n");
        exit(1);
    }

    /*termalizzazione a caldo*/
    rlxd_init(1,123487);
    /*inizializzazione*/
    ranlxd(xx,N);


    /*inizializzo lo stato a sinistra della buca */
    for (t_phys = 0; t_phys < N ; t_phys++)
    {
        xx[t_phys]=shift+(xx[t_phys]-0.5)*2.;
    }

    for ( t_markov = 0; t_markov < 100; t_markov++)
    {
        ranlxd(r3,N);
        ranlxd(r4,N);

        sweep(r3,r4);

        
    }
    /*binning delle posizioni in un istogramma senza decorrelare le posizioni*/
    
    /* 
        Genero i cammini ed eseguo il binning delle posizioni
        infine costruisco uso la funzione di hash per mettere in relazione 
        la posizione della particella al bin di un istogramma
    */

    double Counter;
    for (t_markov=0 ; t_markov < M_sweep; t_markov++)
    {
        ranlxd(r3,N);
        ranlxd(r4,N);
        sweep(r3,r4);
        /*printf("%lf\n", xx[10]);*/
        for (t_phys = 0; t_phys < N ; t_phys++)
        {   
            /*printf("%lf\n", xx[i]);*/
            /*printf("%lf\n", xx[t_phys]);*/
            int Bin = ((int)(0.5+100*(xx[t_phys]+2.5)/(5.)));
            printf("%d\n", Bin);
            /*HISTO[((int)(100*(xx[t_phys]+2.5)/(5.)))]=HISTO[((int)(100*(xx[t_phys]+2.5)/(5.)))]+1.;
            /*printf("%d\n", HISTO[Bin]);*/
            /*printf("%lf\n",HISTO[((int)(100*(xx[t_phys]+2.5)/(5.)))]);*/
            /*fprintf(Histo_Q, "%d\n", Bin);*/
            
            if (Bin<100 && Bin>=0)
            {
                HISTO[Bin] = HISTO[Bin]+1;
                Counter = Counter +1;
                
            }
            
        } 
    }
    
    double a = 1;

    double ABS_x_min = 0;
    /*double x_max = 100;*/
    for ( i = 0; i < 100; i++)
    { 
        /*printf("%lf\n", (double)HISTO[i]/(double)(M_sweep*N));*/
        
        /*calcolo la fdp del primo stato eccitato sottraendo il valore teorico sul reticolo del quadrato 
        della funzione d'onda dello stato di ground (sottraggo |phi_0(x)|^2, 
        dove x è la posizione che corrisponde al centro dle bin )*/

        /*double x =  (x_max+ABS_x_min)/(2*100) + i * ((x_max+ABS_x_min)/100) - ABS_x_min;

        /*printf("%lf\n", x);*/

        /*double X = sqrt(sqrt( (M * W) / 3.14)) * exp( - ((M/2) * W ) * (x*x) );

        /*printf("%lf\n", X);*/
        
        /*printf("%lf\n", X*X);*/
        /*"%lf\n", ((100/4.)*( (double)HISTO[i]/(double)(M_sweep*N))) */
        fprintf(Histo_Q, "%lf ", x_min+i*((x_max-x_min)/100));
        fprintf(Histo_Q, "%lf\n", ((100/5.)*( (double)HISTO[i]/(double)(Counter))));

       

    }
    printf("ratio accettanza");
    printf("%lf\n", (double)Counter/(M_sweep*N));
    
    fclose(Histo_Q);
    
    return 0;


}
