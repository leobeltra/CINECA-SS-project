#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "../include/global.h"
#include "../include/random.h"

#include "../include/sweep_scalar_field.h"
#include "../include/scalar_field_action.h"

int main () {

    /*Create file to record value of the Ising hamitonian*/
    FILE *scalar_action;
    scalar_action=fopen("../../data_analysis/scalar_field_1d/scalar_action.txt","wt");
    if( scalar_action==NULL ) {
        perror("Error in opening file");
        exit(1);
    }

    FILE *field_cnfg;
    field_cnfg=fopen("../../data_analysis/scalar_field_1d/scalar_field_cnfg.txt","wt");
    if( field_cnfg==NULL ) {
        perror("Error in opening file");
        exit(1);
    }
    
    FILE *Histo_Q;
    Histo_Q = fopen("../../data_analysis/scalar_field_1d/hist_state.txt", "wt");
    if (Histo_Q == NULL)
    {
        printf("Problem in opening file.\n");
        exit(1);
    }
    
    FILE *Histo_Q_y;
    Histo_Q_y = fopen("../../data_analysis/scalar_field_1d/hist_state_y.txt", "wt");
    if (Histo_Q_y == NULL)
    {
        printf("Problem in opening file.\n");
        exit(1);
    }

    FILE *Histo_Q_2d;
    Histo_Q_2d = fopen("../../data_analysis/scalar_field_1d/hist_state_2d.txt", "wt");
    if (Histo_Q_2d == NULL)
    {
        printf("Problem in opening file.\n");
        exit(1);
    }
    int time, space;
    double r1[N*N*N];
    double r2[N*N*N];
    double axc;
    double acceptancy;
    int t_markov;
    int rows, cols, i, j;
    int conf=1000;
    int k, y;
    int Counter = 0;
            int Counter_x=0;
        int Counter_y=0;

    double HISTO_2D[101][101];

    int HISTO_x[101];
    int HISTO_y[101];

    for (i = 0; i < 101; i++){   
        HISTO_x[i] = 0;
        HISTO_y[i] = 0;
        for(j=0; j<101; j++){
            HISTO_2D[i][j] = 0;
        }
    }
    

    /*binning variables*/
    int Nbin_max=1000;
    int D_bin=10;

    rlxd_init(1,45321);

    /*initialize field*/
    for (time=0; time<N; time++){

        for (space = 0; space < N; space++){

            for(k=0;k<N; k++){

                phi[time][space][k]=0;

            }
        }
    }

    for(t_markov=0; t_markov<1000; t_markov++){
        /*initialize radom tupled*/
        ranlxd(r1,N*N*N);
        ranlxd(r2,N*N*N);
        
        axc = sweep_scalar_field(r1,r2);

        /*save value of the action at the end of each sweep*/
        fprintf(scalar_action, "%f\n", scalar_field_action());
        /*printf("acceptancy is %lf\n", axc);*/
    }

    double S_mean=0;
    for(t_markov=0; t_markov<conf; t_markov++){
        /*initialize random tuples*/
        ranlxd(r1,N*N*N);
        ranlxd(r2,N*N*N);
        
        /*field configuration*/


        if (t_markov==100)
        {
            for (rows = 5; rows < 6 ; rows++) /*fix at a particular timeslice*/
            {
                for (cols = 0; cols < N; cols++) /*change x*/
                {
                    for(y=0; y<N; y++){

                        fprintf(field_cnfg,"%lf ", phi[rows][cols][y]); /*plot in y*/
                    }
                    
                    fprintf(field_cnfg,"\n");
                }
                
            }
            
        }
        
        axc = sweep_scalar_field(r1,r2);
        acceptancy=acceptancy+axc;
        /*average positions*/
        
        /*for(rows=0; rows<N; rows++){

            for(k=0; k<N; k++){

                for (y=0; y<N;y++){

                    int Bin_y = ((int)(100*(phi[rows][k][y]+5.)/(10.)));
                    
                    if(Bin_y<100 && Bin_y>=0){
                        HISTO_y[Bin_y] = HISTO_y[Bin_y]+1;
                        Counter_y = Counter_y +1;
                    }

                    
                }

   
            }

        }

        for(rows=0; rows<N; rows++){

            for(y=0; y<N; y++){

                for (k=0; k<N;k++){

                    int Bin_x = ((int)(100*(phi[rows][k][y]+5.)/(10.)));
                    
                    if(Bin_x<100 && Bin_x>=0){
                        HISTO_x[Bin_x] = HISTO_x[Bin_x]+1;
                        Counter_x = Counter_x +1;
                    }

                }

   
            }

        }*/
        for(rows=5; rows<6; rows++){

            for(k=0; k<N; k++){

                for (y=0; y<N;y++){

                    int Bin_y = ((int)(100*(phi[rows][k][y]+5.)/(10.)));
                    int Bin_x = ((int)(100*(phi[rows][y][k]+5.)/(10.)));

                    if(Bin_x<100 && Bin_x>=0){
                        HISTO_x[Bin_x] = HISTO_x[Bin_x]+1;
                        Counter_x = Counter_x +1;
                    }
                    
                    if(Bin_y<100 && Bin_y>=0){
                        HISTO_y[Bin_y] = HISTO_y[Bin_y]+1;
                        Counter_y = Counter_y +1;
                    }

                    if (Bin_x<100 && Bin_x>=0){

                        if(Bin_y<100 && Bin_y>=0){                        

                            HISTO_2D[Bin_x][Bin_y] = HISTO_2D[Bin_x][Bin_y]+1;
                            Counter = Counter+1;

                        }
                    }
                }

   
            }

        }


        /*save value of the action at the end of each sweep*/
        fprintf(scalar_action, "%f\n", scalar_field_action());
        S_mean = S_mean + scalar_field_action();
        /*printf("acceptancy is %lf\n", axc);*/
    }
    S_mean = S_mean/conf;
    printf("%lf\n", S_mean);

    acceptancy=acceptancy/conf;
    printf("%lf\n", acceptancy);
    for ( i = 0; i < 100; i++)
    { 
        fprintf(Histo_Q, "%lf\n", ((100./10.)*( (double)HISTO_x[i]/(double)(Counter_x))));
        fprintf(Histo_Q_y, "%lf\n", ((100./10.)*( (double)HISTO_y[i]/(double)(Counter_y))));

        for (j = 0; j < 100; j++)
        {
            fprintf(Histo_Q_2d,"%lf", ((100/10.)*( (double)HISTO_2D[i][j]/(double)(Counter))));
            fprintf(Histo_Q_2d, " ");
        }

        fprintf(Histo_Q_2d,"\n");

    }

    fclose(scalar_action);
    fclose(field_cnfg);
    fclose(Histo_Q);
    fclose(Histo_Q_y);
    fclose(Histo_Q_2d);

    return 0;

}
