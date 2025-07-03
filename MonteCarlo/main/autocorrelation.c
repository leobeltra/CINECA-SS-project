#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/global.h"
#include "../include/sweep.h"
#include "../include/random.h"
#include "../include/azione.h"

int main() {

int t_real;
int t_markov;

int j;
int i;
double y;
int k;
double r3[N];
double r4[N];
int level=1;
int seed=6785;

/*double Y[N][M_sweep];*/

/*malloc matrix with */
double **Y = (double**) malloc(N*sizeof(double*)); 
for(i=0; i<N; i++)
{
	Y[i] = (double*) malloc(M_sweep*sizeof(double));
}

FILE *fb;

/*cold initialization*/
for ( j = 0; j<N; j++) {
    xx[j]=0;
}

/*now we can compute the correlation vector for every sweep */
fb=fopen("../../data_analysis/Gamma_1024.txt","wt");

if( fb==NULL ) {
    perror("Error in opening file");
    exit(1);
}

void rlxd_init(int level,int seed);

/*termalization of the markov chain using 10 times t_markov of termalization */
for ( j = 0; j < 500; j++)
{
    ranlxd(r3,N);
    ranlxd(r4,N);
    sweep(r3,r4);
}

/*after termalization for each sweep C(t) is computed -> it has size = N */
FILE *fp;
fp = fopen("../../data_analysis/correlation_matrix.txt", "wt");
if( fp==NULL ) {
    perror("Error in opening file");
    exit(1);
}

double axc=0;

for ( t_markov = 0; t_markov < M_sweep; t_markov++) {

    ranlxd(r3,N);
    ranlxd(r4,N);
    axc =axc + sweep(r3,r4);

    for (t_real = 0; t_real < N ; t_real++) {
        
        y=0;

        for (k=0; k<N; k++){
            y = y + xx[k]*xx[(k+t_real)%N];
        }

        y=y/N;
        Y[t_real][t_markov] = y; /*correlation matrix*/
        /*fprintf(fp, "%f\n", Y[t_real][t_markov]);*/
    }
}
axc = axc/M_sweep;
printf("l'accettanza è ");
printf("%lf\n", axc);

fclose(fp);

/*compute correlation matrix with row gamma(tm,t)
it will have t=N rows */
/*double** Gamma = (double**)malloc(sizeof(double) *N *M_sweep);*/

double g;
double Gamma[N][M_sweep];
double y_mean;
double Y_mean;
double g_0;
int sweep;
/*computing <y_i(t)>*/
/*compute gamma zero */

/*computing Gamma matrix*/

for (t_real =1; t_real < 2; t_real++) {

    y_mean = 0;
    double Dbin_max = 10000;

    for (t_markov = 0; t_markov < M_sweep; t_markov++) {

        y_mean = y_mean + Y[t_real][t_markov];
    }
    y_mean= y_mean/(M_sweep); 
    /*nota su questo for: quando userò reticoli più grandi dovrò harcodare*/
    /*nota che gamma (T_p = 1 e T_p=63 devono essere uguali: N-i= i )*/
    for (t_markov = 0; t_markov < Dbin_max; t_markov++) {

        g = 0;
        g_0 =0;

        for (sweep = 0; sweep < M_sweep-Dbin_max; sweep++) {

            g  = g + (Y[(t_real)][sweep]*Y[(t_real)][(sweep+t_markov)]);
            g_0 = g_0 +(Y[(t_real)][sweep]*Y[(t_real)][sweep]);
        }
        
        g = g/(M_sweep-Dbin_max)-(y_mean*y_mean);
        g_0 = g_0/(M_sweep-Dbin_max)-(y_mean*y_mean);
        Gamma[t_real][t_markov] = g/g_0;

        fprintf(fb, "%d ", t_markov);
        fprintf(fb, "%lf\n", g/g_0);
    
    }
}

fclose(fb);

FILE *fc;

/*now we can compute the correlation vector for every sweep */
fc=fopen("../../data_analysis/Decorrelation&Binning.txt","wt");

if( fc==NULL ) {
    perror("Error in opening file");
    exit(1);
}

FILE *fd;

/*now we can compute the correlation vector for every sweep */
fd=fopen("../../data_analysis/C_bar.txt","wt");

if( fd==NULL ) {
    perror("Error in opening file");
    exit(1);
}
/*re-compute gamma(t_m,t_phys)/gamma(t_m,0)*/

/*Binning the correlation */
int Dbin = 500;  /*da modificare */

int t_phys;


for (t_phys = 0; t_phys < N; t_phys ++) {
    double C_bar=0;
    double C_tilde=0;
    int Nbin=0;
    /*indice che somma sugli sweep nel bin Nbin-esimo*/
    int j;

    /*scorro le colonne*/
    while (Dbin*Nbin<=M_sweep) {

        for (j=Nbin*Dbin; j < (Nbin+1)*Dbin; j++) {  

            C_tilde = C_tilde + Y[t_phys][j]; 
        }

        C_tilde = C_tilde/Dbin;
        fprintf(fc, "%f\n", C_tilde);
        Nbin = Nbin+1; 
        C_bar = C_bar + C_tilde;
    }

    C_tilde = C_tilde / Nbin;
    fprintf(fd, "%f\n", C_tilde);
}

fclose(fc);
fclose(fd);

return 0;

}
