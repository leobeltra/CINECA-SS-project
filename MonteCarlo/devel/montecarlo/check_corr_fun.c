#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/global.h"
#include "../../include/azione.h"
#include "../../include/dazione.h"
#include "../../include/random.h"
#include "../../include/geometrize_T.h"
#include "../../include/geometrize_V.h"
#include "../../include/geometrize_corr_fun.h"

/**********************************************
*   Check azione e sua variazione
*   Si utilizza la convergenza di serie geometriche
*   con ragione |q|<1.
***********************************************/

int main(){

    rlxd_init(1,78);
    double r1[N];
    ranlxd(r1, N);

    int i;
    double xx_init=1;
    double q=0.5;

    /*Inizializzo il path con una serie geometrica*/

    /*Per farlo scelgo il valore del primo elemento di xx[] e la ragione q*/
    xx[0] = xx_init;
    for (i = 0; i < N; i++)
    {
        xx[i+1]=xx[i]*q;
        /*printf("%lf\n", xx[i]);*/
    }

    /*calcolo il contributo cinetico a S come nel modulo action.c*/
    for ( i = 0; i < N; i++)
    {
        
    }
    /*calcolo il contributo di potenziale a S come nel modulo action.c*/


    /*chiamo le funzione che calcolano la somma parziale di serie geometriche*/
    /*moltiplico il risultato delle somme parziali per le costanti ficiche che compaiono in S*/

    /*somma parziale geometrica per S potenziale*/
    double S_V_geom = geometrize_V(0,q, xx_init, N);
    S_V_geom = S_V_geom*((M*W*W));/*moltiplico per mw^2 che per ora non ho considerato nella somma*/

    /*somma parziale geometrica per il S cinetico*/
    double  S_T_geom = geometrize_T(0,q, xx_init, N)+(xx[N-1]-xx[0])*(xx[N-1]-xx[0])/2.;
    S_T_geom=S_T_geom*(M);/*moltiplico per il termine di massa m/2*/

    /* questo lo devo confrontare con i risultati dati dalla serie geometrica*/
    /* In questo modo siamo sicuri che la routine azione restituisca proprio T+V*/
    double action = azione();   /*questa mi serve solo a verificare che T+V = S*/

    printf(" S_V_geom     S_V       S_T_geom      S_T\n");
    printf("%lf     %lf     %lf     %lf\n", S_V_geom, S_V, S_T_geom, S_T);
    printf(" azione         T+V     T_geom+V_geom\n");
    printf("%lf     %lf     %lf\n", action, S_T+S_V, S_V_geom+S_T_geom);
    
    /*CHECK dS*/
    /*propongo una nuova posizione al tempo fisico fissato per calcolare il dS usando i risultati delle somme geometriche; 
      la variazione viene poi confrontata con la routine dazione.c per verificarne il corretto funzionamento*/

    int j;
    for ( i = 1; i <2; i++)
    {
        double xxnew =xx[i]- 2*D*(r1[i]-0.5);
        double delta_S = -(M/2.)*(q-1)*(q-1)*(xx[i]*xx[i])+((M*W*W)/2.)*(xxnew*xxnew-xx[i]*xx[i]);
        double d_azione = -dazione(xxnew, i);

        printf("%lf\n", (M/2.)*(q-1)*(q-1)*(xxnew*xxnew-xx[i]*xx[i]));
        double S_V_new=0;
        double S_T_new=0;

        xx[i]=xxnew;

        for ( j = 0; j < N; j++)
        {
            S_T_new+=(M/2.0)*(xx[(j+1)%N]-xx[j])*(xx[(j+1)%N]-xx[j]);
        }

        for ( j= 0; j < N; j++)
        {
            S_V_new+=((M*W*W)/2.0)*xx[j]*xx[j];
        }
        printf("%lf\n", S_T_new-S_T);

        printf("%lf %lf %lf\n", delta_S, d_azione, (S_T_new+S_V_new)-(S_T+S_V));
        
    }

    
    


    
    return 0;
}
