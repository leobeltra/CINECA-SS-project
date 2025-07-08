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

__device__
double dH_ising_2D(double s_new, int i, int j, double* ss)
{

    double dH_B = 2*B_field*ss[i*N +j];

    double dH_S = 2*J*ss[i*N +j]*(ss[i*N + ((j+1)%N)]  +ss[i*N+((j+N-1)%N)]   + ss[((i+1)%N)*N+j]  +   ss[((i+N-1)%N)*N+j]);

    double dH_ising=dH_B+dH_S;

    return dH_ising;

}

__global__
void sweep_dH(int i, double* r, double* ss, double* accettanza) {
    int j = blockIdx.x * blockDim.x + threadIdx.x; // calculate the column index
    if (j < N) { // check bounds
        double s_new = -ss[i*N + j]; // flip the spin

        if(exp(-beta*dH_ising_2D(s_new,i,j, ss))>r[(N)*i+j])
        {

            ss[i*N +j]=s_new;
            *accettanza += 1.0; // increment the acceptance count
        } 
        // printf("i=%d j=%d accettanza=%f\n", i, j, *accettanza);
    }    
}

__host__
double sweep_ising_2D(double* r, double* ss, double* accettanza) {

    double s_new;
    int i, j;

    *accettanza = 0.0;

    int blockSize = 32; // numero di thread per blocco
    int gridSize = (N + blockSize - 1) / blockSize; // numero di blocchi per coprire N

    for(i=0; i<N; i++)
    {
        sweep_dH<<<gridSize,blockSize>>>(i, r, ss, accettanza);
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA error after sweep_dH: %s\n", cudaGetErrorString(err));
        }
    }
    cudaDeviceSynchronize(); // Wait for the kernel to finish


    // printf("accettanza: %lf\n", *accettanza);

    *accettanza = *accettanza/(N*N);
    return *accettanza;
}

int main () {

    cudaEvent_t start, stop;
    float elapsedTime;

    // Create events
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Record start event
    cudaEventRecord(start, 0);

double* ss;

cudaMallocManaged(&ss, N * N * sizeof(double));


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
double axc;
int Nbin_max = 1000;        /*bin number*/
int Dbin = 10;               /*bin size*/    
int t_markov;                /*markov time */
int rows, cols;
int Nbin;
double acceptancy=0;

double *accettanza;
cudaMallocManaged(&accettanza, sizeof(double));
cudaMallocManaged(&r, N * N * sizeof(double));

double MAGNETIZATION[Nbin_max]; 

/*double **couplings =(double**) calloc(N);*/

double couplings[2*N][2*N]; 

/*Flat random distribution init*/
rlxd_init(1,263);


/*Initialize matrix*/

/*Chessboard initialization */
for (rows=0; rows<N; rows++){

    for (cols = 0; cols < N; cols++)
    {
        if (rows%2==0)
        {
            ss[rows*N + cols] =1;
        }
        else ss[rows*N + cols] =-1;
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

    axc = sweep_ising_2D(r, ss, accettanza);
    cudaDeviceSynchronize(); // Wait for the kernel to finish

    printf("accettanza: %lf\n", axc);
    fprintf(H_ising, "%f\n", H_Ising_2D(ss)); 

}

printf("Accettanza finale: %lf\n", axc);


while (Nbin*Dbin < M_sweep)
{                       

    for ( t_markov=Nbin*Dbin ; t_markov < (Nbin+1)*Dbin; t_markov++) 
    {
        ranlxd(r,N*N);

        axc=sweep_ising_2D(r, ss, accettanza);
        acceptancy = acceptancy+axc;

        double c=0;
        
        for (rows = 0; rows < N ; rows++)                      
        {   

            for (cols = 0; cols < N; cols++)
            {
                c=c+ss[rows*N + cols];
            }
            
        }

        c=c/(double)(N*N);    

        MAGNETIZATION[Nbin] = MAGNETIZATION[Nbin]+c; 
    }

    Nbin=Nbin+1;     
}

printf("%lf\n", acceptancy/(double)M_sweep);

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

    // Record stop event
    cudaEventRecord(stop, 0);

    // Wait for kernel & events to finish
    cudaEventSynchronize(stop);

    // Compute elapsed time in milliseconds
    cudaEventElapsedTime(&elapsedTime, start, stop);

    printf("Time: %f ms\n", elapsedTime);

return 0;

}
