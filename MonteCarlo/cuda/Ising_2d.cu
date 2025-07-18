#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/global.h"
#include "../include/random.h"

// #include "../include/sweep_ising_2D.h"
// #include "../include/dH_ising_2D.h"
#include "../include/H_ising_2D.h"

__global__
void sweep_ising_2D(double* r, double* ss_new, double* ss_old, double* accettanza) {

    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < N && j < N) { 
        double s = ss_old[i*N + j];  // stato attuale dello spin

        double dH_B = 2*B_field*s;

        double dH_S = 2*J*s*(ss_old[i*N + (j+1)%N]+ss_old[i*N + (j+N-1)%N]+ss_old[((i+1)%N)*N + j]+ss_old[((i+N-1)%N)*N + j]);

        double dH_ising=dH_B+dH_S;

        if (i == 0 && j == 0)
        // printf("s=%f dH=%f prob=%f rand=%f\n", s, dH_ising, exp(-beta * dH_ising), r[i*N + j]);

        if(exp(-beta*dH_ising)>r[(N)*i+j])
        {
            ss_new[i*N+j]=-s;
            atomicAdd(accettanza, 1.0); // Increment the acceptance count atomically
        } else {
            ss_new[i*N + j] = s; 
        }
    }    
}

int main () {

// double **ss = (double**) calloc(N,sizeof(double*));
// for (int i = 0; i < N; i++) {
//     ss[i] = (double*) calloc(N,sizeof(double));
// }
double* ss_new;
cudaMallocManaged(&ss_new, N*N*sizeof(double));
double* ss_old;
cudaMallocManaged(&ss_old, N*N*sizeof(double));
double* r;
cudaMallocManaged(&r, N*N*sizeof(double));
double* accettanza;
cudaMallocManaged(&accettanza, sizeof(double));

/*Create file to record value of the magnetization*/
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

int m;
int k;
// double r[N*N];
double axc;
int Nbin_max = 1000;        /*bin number*/
int Dbin = 10;               /*bin size*/    
int t_markov;                /*markov time */
int rows, cols;
int Nbin;
double acceptancy=0;

double MAGNETIZATION[Nbin_max]; 

/*double **couplings =(double**) calloc(N);*/

double couplings[2*N][2*N]; 

/*Flat random distribution init*/
rlxd_init(1,263);


/*Initialize matrix*/

/*Chessboard initialization */
for (int i = 0; i < N; i++)
  for (int j = 0; j < N; j++)
{
    ss_old[i*N + j] = (drand48() < 0.5) ? -1.0 : 1.0;
    ss_new[i*N + j] = (drand48() < 0.5) ? -1.0 : 1.0;
}   


dim3 blockDim(32, 32);
dim3 gridDim(1, 1);
double* temp;

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

    sweep_ising_2D<<<blockDim, gridDim>>>(r, ss_new, ss_old, accettanza);
    cudaDeviceSynchronize();

    // printf("ss_new[0][0] = %lf\n", ss_new[0*N + 0]);

    *accettanza = *accettanza/(N*N);

    fprintf(H_ising, "%f\n", H_Ising_2D(ss_new));
    printf("ss_old[0] = %f, ss_new[0] = %f\n", ss_old[0], ss_new[0]);

    temp = ss_old;
    ss_old = ss_new;
    ss_new = temp;
    printf("ss_old[0] = %f, ss_new[0] = %f\n", ss_old[0], ss_new[0]);
}



while (Nbin*Dbin < M_sweep)
{                       

    for ( t_markov=Nbin*Dbin ; t_markov < (Nbin+1)*Dbin; t_markov++) 
    {
        ranlxd(r,N*N);

        sweep_ising_2D<<<blockDim, gridDim>>>(r, ss_new, ss_old, accettanza);
        cudaDeviceSynchronize();
        acceptancy = acceptancy+*accettanza;

        double c=0;
        
        for (rows = 0; rows < N ; rows++)                      
        {   

            for (cols = 0; cols < N; cols++)
            {
                c=c+ss_new[rows*N + cols];
            }
            
        }

        c=c/(double)(N*N);    

        MAGNETIZATION[Nbin] = MAGNETIZATION[Nbin]+c; 

        double* temp = ss_old;
        ss_old = ss_new;
        ss_new = temp;    
    }

    Nbin=Nbin+1;     
}

printf("%lf\n", acceptancy/(double)M_sweep);

for (rows = 0; rows < N ; rows++)
{
    for (cols = 0; cols < N; cols++)
    {
        fprintf(ising_state,"%lf ", ss_new[rows*N + cols]);
    }
    fprintf(ising_state,"\n");
}

for (Nbin = 0; Nbin < Nbin_max; Nbin++)
{
    fprintf(magnetization, "%lf\n", MAGNETIZATION[Nbin]/(double)Dbin);
}

double a = sinh(beta*B_field)+(sinh(beta*B_field)*cosh(beta*B_field))/(sqrt(sinh(beta*B_field)*sinh(beta*B_field)+exp(-4*beta*J)));
double b = cosh(beta*B)+sqrt(exp(-4*beta*J)+sinh(beta*B_field)*sinh(beta*B_field));
double magnetization_th = a/b;
printf("theoretical average magnetization is %lf", magnetization_th);

fclose(magnetization);
fclose(H_ising);
fclose(ising_state);

return 0;

}
