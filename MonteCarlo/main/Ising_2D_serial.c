#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../include/global.h"
#include "../include/random.h"

#include "../include/H_ising_2D.h"


double sweep_ising_2D_b(double* r_b, double* ss_w, double* ss_b) {

    double s_new;
    int i, j;
    double accettanza = 0.0;

    for ( i = 0; i < N; i++)
    {
        for ( j = 0; j < (N/2); j++)
        {
            double dH_B = 2*B_field*ss_b[i*N/2+j];

            double dH_S = 2*J*ss_b[i*(N/2)+j]*(ss_w[i*N/2 + ((j+1)%(N/2))]  + ss_w[i*(N/2)+((j+(N/2)-1)%(N/2))]   
                + ss_w[((i+1)%N)*(N/2)+j]  +   ss_w[((i+N-1)%N)*(N/2)+j]);

            double dH_ising=dH_B+dH_S;

            s_new = -ss_b[i*N/2+j];
            if(exp(-beta*dH_ising)>r_b[i*N/2+j])
            {
                ss_b[i*(N/2)+j]=s_new;
                accettanza = accettanza + 1; 
                /* printf("acceptance: %lf\n", accettanza); */
                /* printf( " %lf\n %lf\n, ", ss_b, ss_w); */
            } 
            
        }
        
    }
    accettanza = accettanza/(N*N);
    return accettanza;
    
}

double sweep_ising_2D_w(double* r_w, double* ss_w, double* ss_b) {

    double s_new;
    int i, j;
    double accettanza = 0.0;


    for ( i = 0; i < N; i++)
    {
        for ( j = 0; j < (N/2); j++)
        {
            double dH_B = 2*B_field*ss_w[i*N/2+j];

            double dH_S = 2*J*ss_w[i*(N/2)+j]*(ss_b[i*(N/2) + (j+1)%(N/2)]  +ss_b[i*(N/2)+((j+(N/2)-1)%(N/2))]   + ss_b[((i+1)%N)*(N/2)+j]  +   ss_b[((i+N-1)%N)*(N/2)+j]);

            double dH_ising=dH_B+dH_S;

            s_new = -ss_w[i*(N/2)+j];
            if(exp(-beta*dH_ising)>r_w[i*N/2+j])
            {
                ss_w[i*(N/2)+j]=s_new;
                accettanza = accettanza + 1; 
            } 

        }
        
    }

    accettanza = accettanza/(N*N);
    return accettanza;
    
}


int main () {

        clock_t start, end;
    double cpu_time_used;

    start = clock();    

double *ss_w= (double*) calloc ((N*N)/2, sizeof(double));
double *ss_b= (double*) calloc ((N*N)/2, sizeof(double));
double *r_b = (double*) calloc ((N*N)/2, sizeof(double));
double *r_w = (double*) calloc ((N*N)/2, sizeof(double));

FILE *magnetization;
magnetization=fopen("../../data_analysis/Ising/magnetization_2D.txt","wt");
if( magnetization==NULL ) {
    perror("Error in opening file");
    exit(1);
}

FILE *ising_state;
ising_state=fopen("../../data_analysis/Ising/issing_state.txt","wt");
if( ising_state==NULL ) {
    perror("Error in opening file");
    exit(1);
}

/*Create file to record value of the Ising hamitonian*/
FILE *H_ising;
H_ising=fopen("../../data_analysis/Ising/H_issing_2D.txt","wt");
if( H_ising==NULL ) {
    perror("Error in opening file");
    exit(1);
}

/*variables */
int m;
int k;
double axc;
int Nbin_max = 1000;        /*bin number*/
int Dbin = 10;               /*bin size*/    
int t_markov;                /*markov time */
int rows, cols;
int Nbin;
double acceptancy=0;


double MAGNETIZATION[Nbin_max]; 
 

/*Flat random distribution init*/
rlxd_init(1,263);


/*Initialize matrix*/

/*Chessboard initialization */
for (rows=0; rows<N; rows++){

    for (cols = 0; cols < N/2; cols++)
    {
        ss_w[rows*(N/2)+cols] = 1.0;
        ss_b[rows*(N/2)+cols] = (rand() % 2 == 0) ? 1.0 : -1.0; 
    
    }
}

/*THERMALIZATION OF THE MARKOV CHAIN*/
for (m=0; m<10000; m++)
{
    /* printf("Thermalization step %d\n", m); */
    /*randomize vector r to be used in the Metropolis trials*/
    ranlxd(r_w,N*(N/2));
    ranlxd(r_b,N*(N/2));

    /* ranlxd(r_b,N*N/2); */
    /* ranlxd(r_w,N*N/2); */

    /*if (m==100)
    {
        for (rows = 0; rows < N ; rows++)
        {
            for (cols = 0; cols < N; cols++)
            {
                fprintf(ising_state,"%lf ", ss[rows][cols]);
            }
            fprintf(ising_state,"\n");
        }
          
    }*/

    

    axc = sweep_ising_2D_b(r_b, ss_w, ss_b);
    axc += sweep_ising_2D_w(r_w, ss_w, ss_b);

    fprintf(H_ising, "%f\n", H_ising_2D(ss_b)+H_ising_2D(ss_w)); 

}

printf("Accettanza finale: %lf\n", axc);

/*MARKOV CHAIN MONTECARLO SWEEP*/

    for ( t_markov=0 ; t_markov < M_sweep; t_markov++) 
    {
        ranlxd(r_w,N*(N/2));
        ranlxd(r_b,N*(N/2));

        axc = sweep_ising_2D_b(r_b, ss_w, ss_b);
        
        axc += sweep_ising_2D_w(r_w, ss_w, ss_b);

        acceptancy = acceptancy + axc;

    }



/*MARKOV CHAIN MONTECARLO SWEEP*/

printf("average acceptancy of the MonteCarlo is: %lf\n", acceptancy/(double)M_sweep);

for (rows = 0; rows < N ; rows++)
{
    for (cols = 0; cols < N/2; cols++)
    {
        fprintf(ising_state,"%lf ", ss_b[rows*(N/2) + cols]);
        fprintf(ising_state,"%lf ", ss_w[rows*(N/2) + cols]);
    }
    fprintf(ising_state,"\n");
}

for (Nbin = 0; Nbin < Nbin_max; Nbin++)
{
    fprintf(magnetization, "%lf\n", MAGNETIZATION[Nbin]/(double)Dbin);
}

free(ss_b);
free(ss_w);
free(r_b);
free(r_w);

printf("End of the program\n");

fclose(magnetization);
fclose(H_ising);
fclose(ising_state);

    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Total execution time: %f seconds\n", cpu_time_used);

return 0;

}
