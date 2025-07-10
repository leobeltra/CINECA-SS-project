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
    int id = blockIdx.x * blockDim.x + threadIdx.x;

    int sites_per_row = N/2;  // = 4 per N=8

    int i = id / sites_per_row;  // indice della riga
    int j_in_row = id % sites_per_row;  // posizione nella riga "nera"

    // ora calcoliamo j vero nel lattice 2D
    int j;
    if (i % 2 == 0) {
        j = 2 * j_in_row; // righe pari i=0,2,4,.. iniziano da colonna 1
    } else {
        j = 2 * j_in_row + 1;     // righe dispari i=1,3,5,.. iniziano da colonna 0
    }

    s_new = -ss_b[id];

    double dH_B = 2*B_field*ss_b[id];

    double dH_S = 2*J*ss_b[id]*(ss_w[i*N + ((j+1)%N)]  +ss_w[i*N+((j+N-1)%N)]   + ss_w[((i+1)%N)*N+j]  +   ss_w[((i+N-1)%N)*N+j]);

    double dH_ising=dH_B+dH_S;

    if(exp(-beta*dH_ising)>r[id])
    {
        ss_b[id]=s_new;
        atomicAdd(accettanza, 1.0); // Increment the acceptance count atomically
    } 
}



// white same but with if (i % 2 == 0) j = 2*j_in_row;     // righe pari: bianchi in 0,2,4,6
//                        else            j = 2*j_in_row + 1; // righe dispari: bianchi in 1,3,5,7


__global__
void sweep_ising_2D_w(double* r, double* ss_w, double* ss_b, double* accettanza) {

    double s_new;
    int id = blockIdx.x * blockDim.x + threadIdx.x;

    int sites_per_row = N/2;  // = 4 per N=8

    int i = id / sites_per_row;  // indice della riga
    int j_in_row = id % sites_per_row;  // posizione nella riga "nera"

    // ora calcoliamo j vero nel lattice 2D
    int j;
    if (i % 2 == 0) {
        j = 2 * j_in_row + 1; // righe pari i=0,2,4,.. iniziano da colonna 1
    } else {
        j = 2 * j_in_row;     // righe dispari i=1,3,5,.. iniziano da colonna 0
    }

    s_new = -ss_w[id];

    double dH_B = 2*B_field*ss_w[id];

    double dH_S = 2*J*ss_w[id]*(ss_b[i*N + ((j+1)%N)]  +ss_b[i*N+((j+N-1)%N)]   + ss_b[((i+1)%N)*N+j]  +   ss_b[((i+N-1)%N)*N+j]);

    double dH_ising=dH_B+dH_S;

    if(exp(-beta*dH_ising)>r[id])
    {
        ss_w[id]=s_new;
        *accettanza = *accettanza + 1; // Increment the acceptance count atomically
    } 
}

int main () {

    cudaEvent_t start, stop;
    float elapsedTime;

    // Create events
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Record start event
    cudaEventRecord(start, 0);

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
double* r;
int Nbin_max = 1000;        /*bin number*/
int Dbin = 10;               /*bin size*/    
int t_markov;                /*markov time */
int rows, cols;
int Nbin;
double acceptancy=0;

double *accettanza;
cudaMallocManaged(&accettanza, sizeof(double));
cudaMallocManaged(&r, N * N * sizeof(double) / 2);

double MAGNETIZATION[Nbin_max]; 

/*double **couplings =(double**) calloc(N);*/

double couplings[2*N][2*N]; 

/*Flat random distribution init*/
rlxd_init(1,263);


/*Initialize matrix*/

/*Chessboard initialization */
for (rows=0; rows<N/2; rows++){

    for (cols = 0; cols < N/2; cols++)
    {
        if (rows%2==0)
        {
            ss_w[rows*N/2 + cols] =1;
        }
        else ss_w[rows*N/2 + cols] =-1;
    }
}

for (rows=0; rows<N/2; rows++){

    for (cols = 0; cols < N/2; cols++)
    {
        if (rows%2==0)
        {
            ss_b[rows*N/2 + cols] =1;
        }
        else ss_b[rows*N/2 + cols] =-1;
    }
}

int block_dim = 32;
int grid_dim = (N/2 + block_dim - 1) / block_dim;

//THERMALIZATION OF THE MARKOV CHAIN
for (m=0; m<1000; m++)
{

    /*randomize vector r to be used in the Metropolis trials*/
    ranlxd(r,N*N/2);

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

    sweep_ising_2D_b<<<grid_dim, block_dim>>>(r, ss_w, ss_b, accettanza);
    cudaDeviceSynchronize(); // Wait for the kernel to finish
    sweep_ising_2D_w<<<grid_dim, block_dim>>>(r, ss_w, ss_b, accettanza);
    cudaDeviceSynchronize(); // Wait for the kernel to finish

    *accettanza = *accettanza/(N*N);


    //printf("thermalization acceptancy: %lf\n", *accettanza);
    // fprintf(H_ising, "%f\n", H_Ising_2D(ss)); 

}

//printf("Accettanza finale: %lf\n", *accettanza);

//MARKOV CHAIN MONTECARLO SWEEP
while (Nbin*Dbin < M_sweep)
{                       

    for ( t_markov=Nbin*Dbin ; t_markov < (Nbin+1)*Dbin; t_markov++) 
    {
        ranlxd(r,N*N);

        sweep_ising_2D_b<<<grid_dim, block_dim>>>(r, ss_w, ss_b, accettanza);
        cudaDeviceSynchronize(); // Wait for the kernel to finish
        sweep_ising_2D_w<<<grid_dim, block_dim>>>(r, ss_w, ss_b, accettanza);
        cudaDeviceSynchronize(); // Wait for the kernel to finish

        acceptancy = acceptancy + *accettanza;

        double c=0;
        
        for (rows = 0; rows < N /2; rows++)                      
        {   

            for (cols = 0; cols < N /2; cols++)
            {
                c=c+ss_b[rows*N/2 + cols] + ss_w[rows*N/2 + cols];
            }
            
        }

        c=c/(double)(N*N);    

        MAGNETIZATION[Nbin] = MAGNETIZATION[Nbin]+c; 
    }

    Nbin=Nbin+1;     
}


//MARKOV CHAIN MONTECARLO SWEEP no-decorrelation

printf("average acceptancy of the MonteCarlo is: %lf\n", acceptancy/(double)M_sweep);

// for (rows = 0; rows < N ; rows++)
// {
//     for (cols = 0; cols < N; cols++)
//     {
//         fprintf(ising_state,"%lf ", ss[rows*N + cols]);
//     }
//     fprintf(ising_state,"\n");
// }

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

    // Record stop event
    cudaEventRecord(stop, 0);

    // Wait for kernel & events to finish
    cudaEventSynchronize(stop);

    // Compute elapsed time in milliseconds
    cudaEventElapsedTime(&elapsedTime, start, stop);

    printf("Time: %f ms\n", elapsedTime);

return 0;

}
