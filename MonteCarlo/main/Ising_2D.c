#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/global.h"
#include "../include/random.h"

#include "../include/sweep_ising_2D.h"
#include "../include/dH_ising_2D.h"
#include "../include/H_ising_2D.h"

int main () {

/*Create file to record value of the magnetization*/
FILE *magnetization;
magnetization=fopen("../../data_analysis/Ising/magnetization_2D.txt","wt");
if( magnetization==NULL ) {
    perror("Error in opening file");
    exit(1);
}

FILE *issing_state;
issing_state=fopen("../../data_analysis/Ising/issing_state.txt","wt");
if( issing_state==NULL ) {
    perror("Error in opening file");
    exit(1);
}

/*Create file to record value of the Ising hamitonian*/
FILE *H_issing;
H_issing=fopen("../../data_analysis/Ising/H_issing_2D.txt","wt");
if( H_issing==NULL ) {
    perror("Error in opening file");
    exit(1);
}

int m;
int k;
double r[N*N];
double axc;
int Nbin_max = 1000;        /*bin number*/
int Dbin = 10;               /*bin size*/    
int t_markov;                /*markov time */
int rows, cols;
int Nbin;
double acceptancy=0;

double MAGNETIZATION[Nbin_max]; 

/*double **couplings =(double**) calloc(N);*/

/*Flat random distribution init*/
rlxd_init(1,263);


/*Initialize matrix*/

/*Chessboard initialization */
for (rows=0; rows<N; rows++){

    for (cols = 0; cols < N; cols++)
    {
        ss[rows][cols] = (rand() % 2 == 0) ? 1.0 : -1.0; 
    
    }
}



for (m=0; m<1000; m++)
{

    /*randomize vector r to be used in the Metropolis trials*/
    ranlxd(r,N*N);

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

    axc = sweep_ising_2D(r);

    fprintf(H_issing, "%f\n", H_issing_2D());

}



/*while (Nbin*Dbin < M_sweep)
{*/                       

    for ( t_markov=0 ; t_markov < M_sweep; t_markov++) 
    {
        ranlxd(r,N*N);

        axc=sweep_ising_2D(r);
        acceptancy = acceptancy+axc;

        double c=0;
        
        for (rows = 0; rows < N ; rows++)                      
        {   

            for (cols = 0; cols < N; cols++)
            {
                c=c+ss[rows][cols];
            }
            
        }

        c=c/(double)(N*N);    

        MAGNETIZATION[Nbin] = MAGNETIZATION[Nbin]+c; 
    }

    /*Nbin=Nbin+1;     
}*/

printf("acceptancy is: %lf\n", acceptancy/(double)M_sweep);

for (rows = 0; rows < N ; rows++)
{
    for (cols = 0; cols < N; cols++)
    {
        fprintf(issing_state,"%lf ", ss[rows][cols]);
    }
    fprintf(issing_state,"\n");
}

for (Nbin = 0; Nbin < Nbin_max; Nbin++)
{
    fprintf(magnetization, "%lf\n", MAGNETIZATION[Nbin]/(double)Dbin);
}

double a = sinh(beta*B_field)+(sinh(beta*B_field)*cosh(beta*B_field))/(sqrt(sinh(beta*B_field)*sinh(beta*B_field)+exp(-4*beta*J)));
double b = cosh(beta*B)+sqrt(exp(-4*beta*J)+sinh(beta*B_field)*sinh(beta*B_field));
double magnetization_th = a/b;
/*printf("theoretical average magnetization is %lf", magnetization_th);*/

fclose(magnetization);
fclose(H_issing);
fclose(issing_state);

return 0;

}
