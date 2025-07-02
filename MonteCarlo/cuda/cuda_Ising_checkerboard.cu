#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../include/global.h"
#include "../include/random.h"

// #include "../include/sweep_ising_2D.h"
// #include "../include/dH_ising_2D.h"
#include "../include/H_ising_2D.h"

#define BLOCK_SIZE 16

__global__
void sweep_ising_2D_b(double* r, double* ss, double* accettanza) {

    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N || j >= N) return;

    if ((i + j) % 2 == 1) return;

    // // Dimensioni del blocco + margini
    // int tx = threadIdx.x;
    // int ty = threadIdx.y;
    // __shared__ double tile[BLOCK_SIZE+2][BLOCK_SIZE+2];

    // // Carica spin da global a shared (con margini)
    // int global_i = i;
    // int global_j = j;
    
    // // Posizione nel tile (offset di +1 per margini)
    // int local_i = ty + 1;
    // int local_j = tx + 1;

    // // Copia spin centrale
    // tile[local_i][local_j] = ss[global_i * N + global_j];

    // // Carica margini (solo i thread sul bordo lo fanno)
    // if (tx == 0)   tile[local_i][0]        = ss[global_i * N + (global_j + N - 1) % N];
    // if (tx == blockDim.x - 1) tile[local_i][BLOCK_SIZE+1] = ss[global_i * N + (global_j + 1) % N];

    // if (ty == 0)   tile[0][local_j]        = ss[((global_i + N - 1) % N) * N + global_j];
    // if (ty == blockDim.y - 1) tile[BLOCK_SIZE+1][local_j] = ss[((global_i + 1) % N) * N + global_j];

    // // Barriera: aspetta che tutti i thread abbiano copiato
    // __syncthreads();

    // // Leggi spin e vicini da shared memory
    // double s = tile[local_i][local_j];

    // double dH_B = 2 * B_field * s;
    // double dH_S = 2 * J * s * (
    //     tile[local_i][local_j + 1] +
    //     tile[local_i][local_j - 1] +
    //     tile[local_i + 1][local_j] +
    //     tile[local_i - 1][local_j]);

    // double dH_ising = dH_B + dH_S;

    double s = ss[i*N + j];  // stato attuale dello spin

    double dH_B = 2*B_field*s;

    double dH_S = 2*J*s*(ss[i*N + (j+1)%N]+ss[i*N + (j+N-1)%N]+ss[((i+1)%N)*N + j]+ss[((i+N-1)%N)*N + j]);

    double dH_ising=dH_B+dH_S;
    // printf("s=%f dH=%f prob=%f rand=%f\n", s, dH_ising, exp(-beta * dH_ising), r[i*N + j]);

    if(exp(-beta*dH_ising)>r[(N)*i+j])
    {
        ss[i*N+j]=-s;
        atomicAdd(accettanza, 1.0);; // Increment the acceptance count atomically
    }    
}

__global__
void sweep_ising_2D_w(double* r, double* ss, double* accettanza) {

    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N || j >= N) return;

    if ((i + j) % 2 == 0) return;

    double s = ss[i*N + j];  // stato attuale dello spin

    double dH_B = 2*B_field*s;

    double dH_S = 2*J*s*(ss[i*N + (j+1)%N]+ss[i*N + (j+N-1)%N]+ss[((i+1)%N)*N + j]+ss[((i+N-1)%N)*N + j]);

    double dH_ising=dH_B+dH_S;

    // printf("s=%f dH=%f prob=%f rand=%f\n", s, dH_ising, exp(-beta * dH_ising), r[i*N + j]);

    if(exp(-beta*dH_ising)>r[(N)*i+j])
    {
        ss[i*N+j]=-s;
        atomicAdd(accettanza, 1.0);; // Increment the acceptance count atomically
    }
}


int main () {

// double **ss = (double**) calloc(N,sizeof(double*));
// for (int i = 0; i < N; i++) {
//     ss[i] = (double*) calloc(N,sizeof(double));
// }
struct timespec start, stop;
double elapsed_ms;

// Avvia il cronometro
clock_gettime(CLOCK_MONOTONIC, &start);

double* ss;
cudaMallocManaged(&ss, N*N*sizeof(double));
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

int device;
cudaGetDevice(&device);

/*double **couplings =(double**) calloc(N);*/

double couplings[2*N][2*N]; 

/*Flat random distribution init*/
rlxd_init(1,263);


/*Initialize matrix*/

/*Chessboard initialization */
for (int i = 0; i < N; i++)
  for (int j = 0; j < N; j++)
{
    ss[i*N+j] = (rand() % 2 == 0) ? 1.0 : -1.0; 
}   


dim3 blockDim(BLOCK_SIZE, BLOCK_SIZE);
dim3 gridDim((N + BLOCK_SIZE - 1) / BLOCK_SIZE, (N + BLOCK_SIZE - 1) / BLOCK_SIZE);

for (m=0; m<1000; m++)
{

    /*randomize vector r to be used in the Metropolis trials*/
    ranlxd(r,N*N);

    cudaMemset(accettanza, 0, sizeof(double));// Reset acceptance count for each sweep

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

    // Prefetch + sincronizza
    cudaMemPrefetchAsync(r, N*N*sizeof(double), device, 0);
    cudaMemPrefetchAsync(ss, N*N*sizeof(double), device, 0);
    cudaMemPrefetchAsync(accettanza, sizeof(double), device, 0);
    // cudaDeviceSynchronize();

    sweep_ising_2D_b<<<blockDim, gridDim>>>(r, ss, accettanza);
    cudaDeviceSynchronize();
    sweep_ising_2D_w<<<blockDim, gridDim>>>(r, ss, accettanza);
    cudaDeviceSynchronize();

    // printf("ss_new[0][0] = %lf\n", ss_new[0*N + 0]);

    // Torna al CPU se accedi da host
    cudaMemPrefetchAsync(ss, N*N*sizeof(double), cudaCpuDeviceId, 0);
    // cudaDeviceSynchronize();

    fprintf(H_ising, "%f\n", H_issing_2D(ss));
    /* printf("ss_old[0] = %f, ss_new[0] = %f\n", ss);

    // temp = ss_old;
    // ss_old = ss_new;
    // ss_new = temp;
    printf("ss_old[0] = %f, ss_new[0] = %f\n", ss_old[0], ss_new[0]); */
}

*accettanza = *accettanza/(N*N);

printf("Acceptance rate after equilibration: %lf\n", *accettanza);

for ( t_markov=0 ; t_markov < M_sweep; t_markov++) {
    ranlxd(r,N*N);

    cudaMemset(accettanza, 0, sizeof(double)); // Reset acceptance count for each sweep

    // Prefetch + sincronizza
    cudaMemPrefetchAsync(r, N*N*sizeof(double), device, 0);
    cudaMemPrefetchAsync(ss, N*N*sizeof(double), device, 0);
    cudaMemPrefetchAsync(accettanza, sizeof(double), device, 0);
    // cudaDeviceSynchronize();

    sweep_ising_2D_b<<<blockDim, gridDim>>>(r, ss, accettanza);
    cudaDeviceSynchronize();
    sweep_ising_2D_w<<<blockDim, gridDim>>>(r, ss, accettanza);
    cudaDeviceSynchronize();
    *accettanza = *accettanza/(N*N);

    //     double c=0;
        
    //     for (rows = 0; rows < N ; rows++)                      
    //     {   

    //         for (cols = 0; cols < N; cols++)
    //         {
    //             c=c+ss[rows*N + cols];
    //         }
            
    //     }

    //     c=c/(double)(N*N);    

    //     MAGNETIZATION[Nbin] = MAGNETIZATION[Nbin]+c; 

    //     /*
    //     double* temp = ss_old;
    //     ss_old = ss_new;
    //     ss_new = temp;
    //     */    
    // }

    // Torna al CPU se accedi da host
    // cudaMemPrefetchAsync(ss, N*N*sizeof(double), cudaCpuDeviceId, 0);
    // cudaDeviceSynchronize();
    
    // Nbin=Nbin+1;     
    acceptancy = acceptancy + *accettanza;

}

printf("average acceptancy of the MonteCarlo is: %lf\n", acceptancy/(double)M_sweep);

for (rows = 0; rows < N ; rows++)
{
    for (cols = 0; cols < N; cols++)
    {
        fprintf(ising_state,"%lf ", ss[rows*N + cols]);
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
printf("theoretical average magnetization is %lf\n", magnetization_th);


fclose(magnetization);
fclose(H_ising);
fclose(ising_state);

// Ferma il cronometro
clock_gettime(CLOCK_MONOTONIC, &stop);

// Calcola il tempo in millisecondi
elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
elapsed_ms += (stop.tv_nsec - start.tv_nsec) / 1000000.0;

printf("Tempo GPU (C): %f ms\n", elapsed_ms);

return 0;

}
