#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/global.h"
#include "../include/random.h"

#include "../include/sweep_ising.h"
#include "../include/dH_ising.h"
#include "../include/H_Ising.h"

int main () {

/*Create file to record value of the Ising hamitonian*/
FILE *H_ising;
H_ising=fopen("../../data_analysis/Ising/H_ising.txt","wt");
if( H_ising==NULL ) {
    perror("Error in opening file");
    exit(1);
}

/*Create file to record value of the magnetization*/
FILE *magnetization;
magnetization=fopen("../../data_analysis/Ising/magnetization.txt","wt");
if( magnetization==NULL ) {
    perror("Error in opening file");
    exit(1);
}

int m;
int k;
double r[N];
double axc;
int Nbin_max = 100;        /*bin number*/
int Dbin = 10000;          /*bin size*/
int t_markov;              /*markov time */
int lattice_site;
int Nbin;
double MAGNETIZATION[Nbin_max]; 

/*Flat random distribution init*/
rlxd_init(1,4672);

/*Initialization of the Ising lattice with s[i]=1 for each site*/
for (k=0; k<N; k++){
    s[k] =-1;
}

for (m=0; m<M_sweep; m++){

    /*randomize vector r to be used in the Metropolis trial*/
    ranlxd(r,N);

    axc = sweep_ising(r);
    
    fprintf(H_ising, "%f\n", H_Ising());

    /*printf("acceptancy is %lf\n", axc);*/
}


while (Nbin*Dbin < M_sweep)
{                       


    for ( t_markov=Nbin*Dbin ; t_markov < (Nbin+1)*Dbin; t_markov++) 
    {
        ranlxd(r,N);

        sweep_ising(r);

        double c=0;

        for (lattice_site = 0; lattice_site < N ; lattice_site++)                      
        {   
            int k;

            c=c+s[(lattice_site)];
        }

        c=c/(double)N;    

        MAGNETIZATION[Nbin] = MAGNETIZATION[Nbin]+c; 
    }

    Nbin=Nbin+1;     
}

for (Nbin = 0; Nbin < Nbin_max; Nbin++)
{
    fprintf(magnetization, "%lf\n", MAGNETIZATION[Nbin]/(double)Dbin);
}


double a = sinh(beta*B_field)+(sinh(beta*B_field)*cosh(beta*B_field))/(sqrt(sinh(beta*B_field)*sinh(beta*B_field)+exp(-4*beta*J)));
double b = cosh(beta*B)+sqrt(exp(-4*beta*J)+sinh(beta*B_field)*sinh(beta*B_field));
double magnetization_th = a/b;
printf("theoretical average magnetization is %lf", magnetization_th);

return 0;

}
