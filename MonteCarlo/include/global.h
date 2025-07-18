
/*******************************************************************************
*
* File global.h
*
* Global parameters and arrays
*
*******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

/*parameters for QFT free and quadratic theory*/
#define M 1
#define W 1

/*lattice spacing*/
#define SPACING 1

/*radius of flat distribution*/
#define D 1

#define lambda 0.01
#define E 1

/*scalar field QFT coupling -  quartic iteraction QM coupling*/
#define G 0.01

#define K 1

/*Other couplings i do not remember*/
#define A 1
#define B 0.25
#define P 0

/*Ising model parameters*/
#define J 0.1
#define B_field 0.0
#define beta 4.5

/*lattice size*/
#define N 256

/*Monte Carlo steps - Markovian time*/
#define M_sweep 10000

#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif

/*lattices*/

/* EXTERN double xx[N]; */


/* EXTERN double yy[N]; */


/* EXTERN double s[N]; */


/* EXTERN double ss[N][N]; */

/* EXTERN double phi[N][N][N]; */


/*EXTERN double sg[N][N][N];*/

#undef EXTERN

#endif

