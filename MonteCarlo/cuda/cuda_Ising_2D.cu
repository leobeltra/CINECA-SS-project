#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/global.h"
#include "../include/random.h"

#include "../include/H_ising_2D.h"

// __global__ 
// void print_info() {
//     printf("This is a CUDA kernel function.\n");
// }

// __device__
// double dH_ising_2D_b(int i, int j, double* ss_w, double* ss_b)
// {

//     double dH_B = 2*B_field*ss_b[i*N +j];

//     double dH_S = 2*J*ss_b[i*N +j]*(ss_w[i*N + ((j+1)%N)]  +ss_w[i*N+((j+N-1)%N)]   + ss_w[((i+1)%N)*N+j]  +   ss_w[((i+N-1)%N)*N+j]);

//     double dH_ising=dH_B+dH_S;

//     return dH_ising;

// }

__global__
void sweep_ising_2D_b(double* r, double* ss_w, double* ss_b, double* accettanza) {

    double s_new;
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < N && j < N/2) { 

        s_new = ss_b[i*(N/2) + j];  // ss_b is a 1D array of size N*N/2, so we access it with i*(N/2) + j
        // printf("id: %d, i: %d, j: %d, s_new: %lf\n", id, i, j, s_new);
    
        ////////////////////////////////////////////////////////////////////
        double dH_B = 2*B_field*s_new;

        double dH_S = 2*J*s_new*(ss_w[i*(N/2) + ((j+1)%(N/2))]  +ss_w[i*(N/2)+((j+N/2-1)%(N/2))]   + ss_w[((i+1)%N)*(N/2)+j]  +   ss_w[((i+N-1)%N)*(N/2)+j]);

        double dH_ising=dH_B+dH_S;

        if(exp(-beta*dH_ising)>r[i*(N/2) + j])
        {
            ss_b[i*(N/2) +j]=-s_new;
            atomicAdd(accettanza, 1.0);  // Increment the acceptance count atomically
        } 
    }
}


// white same but with if (i % 2 == 0) j = 2*j_in_row;     // righe pari: bianchi in 0,2,4,6
//                        else            j = 2*j_in_row + 1; // righe dispari: bianchi in 1,3,5,7


__global__
void sweep_ising_2D_w(double* r, double* ss_w, double* ss_b, double* accettanza) 
{
    double s_new;
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < N && j < N/2) { 

        s_new = ss_w[i*(N/2) + j];  // ss_b is a 1D array of size N*N/2, so we access it with i*(N/2) + j
        // printf("id: %d, i: %d, j: %d, s_new: %lf\n", id, i, j, s_new);
    
        ////////////////////////////////////////////////////////////////////
        double dH_B = 2*B_field*s_new;

        double dH_S = 2*J*s_new*(ss_b[i*(N/2) + ((j+1)%(N/2))]  +ss_b[i*(N/2)+((j+N/2-1)%(N/2))]   + ss_b[((i+1)%N)*(N/2)+j]  +   ss_b[((i+N-1)%N)*(N/2)+j]);

        double dH_ising=dH_B+dH_S;

        if(exp(-beta*dH_ising)>r[i*(N/2) + j])
        {
            ss_w[i*(N/2) +j]=-s_new;
            atomicAdd(accettanza, 1.0);  // Increment the acceptance count atomically
        } 
    }
}

int main () 
{

double* ss_b, *ss_w;

cudaMallocManaged(&ss_w, N * N * sizeof(double) / 2);
cudaMallocManaged(&ss_b, N * N * sizeof(double) / 2);

// double *ss = (double**) calloc(N,sizeof(double*));
// for (int i = 0; i < N; i++) {
//     ss[i] = (double*) calloc(N,sizeof(double));
// }


    // int block_dim = 4;
    // // int grid_dim = (N + block_dim - 1) / block_dim;


    // print_info<<<1, block_dim>>>();
    // cudaDeviceSynchronize();

    

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
double* r_b, *r_w;
int Nbin_max = 1000;        /*bin number*/
int Dbin = 10;               /*bin size*/    
int t_markov;                /*markov time */
int rows, cols;
int Nbin;
double acceptancy=0;

double *accettanza;
cudaMallocManaged(&accettanza, sizeof(double));
cudaMallocManaged(&r_w, N * N * sizeof(double) / 2);
cudaMallocManaged(&r_b, N * N * sizeof(double) / 2);

double MAGNETIZATION[Nbin_max]; 

/*double **couplings =(double**) calloc(N);*/

double couplings[2*N][2*N]; 

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

dim3 block_dim(32, 32);
dim3 grid_dim((N/2 + block_dim.x - 1)/block_dim.x,
              (N   + block_dim.y - 1)/block_dim.y);
//THERMALIZATION OF THE MARKOV CHAIN
for (m=0; m<1000; m++)
{

    /*randomize vector r to be used in the Metropolis trials*/
    ranlxd(r_w,N*(N/2));
    ranlxd(r_b,N*(N/2));

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
    *accettanza = 0.0; 

    sweep_ising_2D_b<<<grid_dim, block_dim>>>(r_b, ss_b, ss_w, accettanza);
    cudaDeviceSynchronize(); // Wait for the kernel to finish

    sweep_ising_2D_w<<<grid_dim, block_dim>>>(r_w, ss_b, ss_w, accettanza);
    cudaDeviceSynchronize(); // Wait for the kernel to finish

    *accettanza = *accettanza/(N*N);


    //printf("thermalization acceptancy: %lf\n", *accettanza);
    fprintf(H_ising, "%f\n", H_Ising_2D(ss_b)+H_Ising_2D(ss_w)); 
}

printf("Accettanza finale: %lf\n", *accettanza);

//MARKOV CHAIN MONTECARLO SWEEP
// while (Nbin*Dbin < M_sweep)
// {                       

     for ( t_markov=0 ; t_markov < M_sweep; t_markov++) 
     {
        ranlxd(r_w,N*N/2);
        ranlxd(r_b,N*N/2);

        *accettanza = 0.0;

        sweep_ising_2D_b<<<grid_dim, block_dim>>>(r_b, ss_b, ss_w, accettanza);
        cudaDeviceSynchronize(); // Wait for the kernel to finish
        sweep_ising_2D_w<<<grid_dim, block_dim>>>(r_w, ss_b, ss_w, accettanza);
        cudaDeviceSynchronize(); // Wait for the kernel to finish

        *accettanza = *accettanza/(N*N);

        acceptancy = acceptancy + *accettanza;

     }



//MARKOV CHAIN MONTECARLO SWEEP

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

//double a = sinh(beta*B_field)+(sinh(beta*B_field)*cosh(beta*B_field))/(sqrt(sinh(beta*B_field)*sinh(beta*B_field)+exp(-4*beta*J)));
//double b = cosh(beta*B)+sqrt(exp(-4*beta*J)+sinh(beta*B_field)*sinh(beta*B_field));
//double magnetization_th = a/b;
//printf("theoretical average magnetization is %lf\n", magnetization_th);

fclose(magnetization);
fclose(H_ising);
fclose(ising_state);

cudaFree(ss_b);
cudaFree(ss_w);
cudaFree(r_b);
cudaFree(r_w);
cudaFree(accettanza);

return 0;

}
