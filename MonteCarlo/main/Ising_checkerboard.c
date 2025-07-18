#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/global.h"
#include "../include/random.h"
#include "../include/H_ising_2D.h"

/*black grid kernel*/


double sweep_ising_2D_b(double** r, double** ss) {

    double s_new;
    int i, j;
    int j_init =0;
    double accettanza = 0.0;

    for ( i = 0; i < N; i++)
    {
        for ( j= j_init; j < N; j += 2)
        {
            double dH_B = 2*B_field*ss[i][j];

            double dH_S = 2*J*ss[i][j]*(ss[i][(j+1)%N] +ss[i][(j+N-1)%N]  + ss[(i+1)%N][j]  +   ss[(i+N-1)%N][j]);

            double dH_ising=dH_B+dH_S;

            s_new = -ss[i][j];

            if(exp(-beta*dH_ising)>r[i][j])
            {
                ss[i][j]=s_new;
                accettanza = accettanza + 1; 
                /* printf("acceptance: %lf\n", accettanza); */
                /* printf( " %lf\n %lf\n, ", ss_b, ss_w); */
            } 

        }

        j_init = (j_init+1)%2;
    }

    accettanza = accettanza/(N*N);
    return accettanza;
    
}



/*white grid kernel*/  
double sweep_ising_2D_w(double** r, double** ss) {

    double s_new;
    int i, j;
    int j_init =1;
    double accettanza = 0.0;

    for ( i = 0; i < N; i++)
    {
        for ( j= j_init; j < N; j += 2)
        {
            double dH_B = 2*B_field*ss[i][j];

            double dH_S = 2*J*ss[i][j]*(ss[i][(j+1)%N] +ss[i][(j+N-1)%N]  + ss[(i+1)%N][j]  +   ss[(i+N-1)%N][j]);

            double dH_ising=dH_B+dH_S;

            s_new = -ss[i][j];

            if(exp(-beta*dH_ising)>r[i][j])
            {
                ss[i][j]=s_new;
                accettanza = accettanza + 1; 
                /* printf("acceptance: %lf\n", accettanza); */
                /* printf( " %lf\n %lf\n, ", ss_b, ss_w); */
            } 

        }

        j_init = (j_init+1)%2;
    }

    accettanza = accettanza/(N*N);
    return accettanza;
    
}

int main () {


/*Flat random distribution init*/
rlxd_init(1,263);


double **ss= (double**) calloc (N, sizeof(double*));
double **r = (double**) calloc (N, sizeof(double*));

for (int i = 0; i < N; i++)
{
    ss[i] = (double*) calloc(N,sizeof(double));
    r[i] = (double*) calloc(N,sizeof(double));
}

FILE *magnetization;
magnetization=fopen("../../data_analysis/Ising/magnetization_2D.txt","wt");
if( magnetization==NULL ) {
    perror("Error in opening file");
    exit(1);
}

FILE *ising_state;
ising_state=fopen("../../data_analysis/Ising/ising_state.txt","wt");
if( ising_state==NULL ) {
    perror("Error in opening file");
    exit(1);
}

/*Create file to record value of the Ising hamitonian*/
FILE *H_ising;
H_ising=fopen("../../data_analysis/Ising/H_ising_2D.txt","wt");
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
 


/*Initialize matrix*/

/*Chessboard initialization */
for (rows=0; rows<N; rows++){

    for (cols = 0; cols < N; cols++)
    {
        ss[rows][cols] = (rand() % 2 == 0) ? 1.0 : -1.0; 
        }
}

printf("initial energy is %f\n", H_ising_2D(ss)); 


/*THERMALIZATION OF THE MARKOV CHAIN*/
for (m=0; m<1000; m++)
{
    /* printf("Thermalization step %d\n", m); */
    /*randomize vector r to be used in the Metropolis trials*/
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            ranlxd(&r[i][j], 1);
        }
    }

    axc = sweep_ising_2D_b(r, ss);
    axc += sweep_ising_2D_w(r, ss);

    fprintf(H_ising, "%lf\n", H_ising_2D(ss));


}

printf("Accettanza finale: %lf\n", axc);

/*MARKOV CHAIN MONTECARLO SWEEP*/

    for ( t_markov=0 ; t_markov < M_sweep; t_markov++) 
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                ranlxd(&r[i][j], 1);
            }
        }


        axc = sweep_ising_2D_b(r, ss);
        axc += sweep_ising_2D_w(r, ss);

        acceptancy = acceptancy + axc;

    }

printf("average acceptancy of the MonteCarlo is: %lf\n", acceptancy/(double)M_sweep);

/*print Ising grid*/

for (rows = 0; rows < N ; rows++)
{
    for (cols = 0; cols < N; cols++)
    {
        fprintf(ising_state,"%lf ", ss[rows][cols]);
    }

    fprintf(ising_state,"\n");
}

/*print average magnetization*/

for (Nbin = 0; Nbin < Nbin_max; Nbin++)
{
    fprintf(magnetization, "%lf\n", MAGNETIZATION[Nbin]/(double)Dbin);
}

free(ss);
free(r);

printf("End of the program\n");

fclose(magnetization);
fclose(H_ising);
fclose(ising_state);

return 0;

}
