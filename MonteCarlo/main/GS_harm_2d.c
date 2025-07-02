/*
############################  @ GS_HARM_2D @   ################################
#
# Federico De Matteis, Lab. di Fisica Computazionale 2022/2023
# Stato di vuoto dell'oscillatore armonico bidimensionale.
# 
##################################################################################
*/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/global.h"
#include "../include/sweep.h"
#include "../include/sweep_2d.h"
#include "../include/random.h"
#include "../include/dazione.h"
#include "../include/dazione_2d.h"
#include "../include/dazione_2d_an.h"
#include "../include/azione.h"

int main()
{
    int HISTO_x[101];
    int HISTO_y[101];
    int HISTO_2D[101][101];
    int i, j;

    for (i = 0; i < 101; i++)
    {   
        HISTO_x[i] = 0;
        HISTO_y[i] = 0;
        for (j = 0; j < 101; j++)
        {
            HISTO_2D[i][j] = 0;
        }
    }

    int t_phys;
    int t_markov;
    double r1[3*N];
    double r2[N];
    double r3[N];

    int Counter=0;

    FILE *Histo_Q_x;
    Histo_Q_x = fopen("../../data_analysis/Histo_Q_x.txt", "wt");
    if (Histo_Q_x == NULL)
    {
        printf("Problem in opening file.\n");
        exit(1);
    }
     FILE *Histo_Q_y;
    Histo_Q_y = fopen("../../data_analysis/Histo_Q_y.txt", "wt");
    if (Histo_Q_y == NULL)
    {
        printf("Problem in opening file.\n");
        exit(1);
    }

    FILE *Histo_Q_2d;
    Histo_Q_2d = fopen("../../data_analysis/Histo_Q_2d.txt", "wt");
    if (Histo_Q_2d == NULL)
    {
        printf("Problem in opening file.\n");
        exit(1);
    }
   
    /*termalizzazione a caldo*/
    rlxd_init(1,1834);
    ranlxd(xx,N);
    ranlxd(yy,N);    
    
    for ( t_markov = 0; t_markov < 5000; t_markov++)
    {
        ranlxd(r1,3*N);
        /*ranlxd(r2,N);
        ranlxd(r3,N);*/
        sweep_2d(r1);
    }
    
    /* 
        Genero i cammini ed eseguo il binning delle posizioni in 2 dimensioni
        infine costruisco una funzione di hash per mettere in relazione 
        la posizione della particella al bin di un istogramma in 2 dimensioni.
    */

    for (t_markov=0 ; t_markov < M_sweep; t_markov++)
    {
        ranlxd(r1,3*N);
        /*ranlxd(r2,N);
        ranlxd(r3,N);*/

        sweep_2d(r1);
 
        for (t_phys = 0; t_phys < N ; t_phys++)
        {   
            /*printf("%lf\n", xx[t_phys]);*/
            int Bin_x = ((int)(100*(xx[t_phys]+2.5)/(5.)));
            int Bin_y = ((int)(100*(yy[t_phys]+2.5)/(5.)));
            
            /*printf("%d\n", Bin_x);
            printf("%d\n", Bin_y);

            /*printf("%d\n", xx[t_phys]);
            printf("%d\n", yy[t_phys]);*/

            if (Bin_x<100 && Bin_x>=0)
            {
                if (Bin_y<100 && Bin_y>=0)
                {
                    HISTO_x[Bin_x] = HISTO_x[Bin_x]+1;
                    HISTO_y[Bin_y] = HISTO_y[Bin_y]+1;

                    Counter = Counter +1;

                    HISTO_2D[Bin_x][Bin_y] = HISTO_2D[Bin_x][Bin_y] + 1; 
                }  
            }
        } 
    }
    
    for ( i = 0; i < 100; i++)
    { 
        fprintf(Histo_Q_x, "%lf\n", ((100/5.)*( (double)HISTO_x[i]/(double)(M_sweep*N))));
        fprintf(Histo_Q_y, "%lf\n", ((100/5.)*( (double)HISTO_y[i]/(double)(M_sweep*N))));

        for (j = 0; j < 100; j++)
        {
            fprintf(Histo_Q_2d,"%lf", ((100/5.)*( (double)HISTO_2D[i][j]/(double)(Counter))));
            fprintf(Histo_Q_2d, " ");
        }

        fprintf(Histo_Q_2d,"\n");
    }
    
    fclose(Histo_Q_2d);
    fclose(Histo_Q_y);
    fclose(Histo_Q_x);
    
    return 0;


}
