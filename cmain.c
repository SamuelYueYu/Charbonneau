#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

extern void init_random(int seed,int seed2);
extern double Xrandom(void);
extern bool cRRG(int conn, int N, int **V);
extern void find_neighbors(int conn, int N, int **V, int *quantiNN, int *quantiNNN, int **NN, int **NNN);

double initialization(int conn, int N, int **V, int *quantiNN, int *quantiNNN, int **NN, int **NNN, double K1, double K2, int *S){
    int i,s1;
    double E;
    
    for(i=0;i<N;i++)
    {
        if(Xrandom()<0.5) S[i] = 1;
        else S[i] = -1;
    }
    E = 0.;
    for(i=0;i<N;i++)
    {
        // the number of neighbors should be revised according to varying c
        // to be revised: 5 -> conn
        for(s1=0;s1<5;s1++) E-= S[i]*S[V[i][s1]];   // nearest neighbors
        for(s1=0;s1<quantiNN[i];s1++) E+= K1*S[i]*S[NN[i][s1]]; // next-nearest neighbors
        for(s1=0;s1<quantiNNN[i];s1++) E+= K2*S[i]*S[NNN[i][s1]];   // next-next-nearest neighbors
    }
    return E/2.;
}

void MCstep(int conn, int N, int **V, int *quantiNN, int *quantiNNN, int **NN, int **NNN, double K1, double K2, double T, int *S, double *E, int *S0, double *C){
    int i,Sflip,sum1,sum2,sum3,s;
    bool accetta;
    double DeltaE,p;
    
    i = (int) (N*Xrandom());
    Sflip = -S[i];
    sum1 = 0; sum2 = 0; sum3 = 0;
    for(s=0;s<conn;s++) sum1+= S[V[i][s]];
    for(s=0;s<quantiNN[i];s++) sum2+= S[NN[i][s]];
    for(s=0;s<quantiNNN[i];s++) sum3+= S[NNN[i][s]];
    DeltaE = (Sflip - S[i])*(-sum1 + K1*sum2 + K2*sum3);
    if(DeltaE<0) accetta = true;
    else{
        p = Xrandom();
        if(p<exp(-DeltaE/T)) accetta = true;
        else accetta = false;
    }
    if(accetta){
        *C+= (Sflip - S[i])*S0[i];
        S[i] = Sflip;
        *E+= DeltaE;
    }
}

int main(int argc, char **argv)
{
    int N,i,**V,*quantiNN,*quantiNNN,**NN,**NNN,*S,j,T_w,nmc,*S0,quanti_punti,index,*time,NRRG,ND,cnt,r,conn,quanti_out;
    double T,K1,K2,kappa,E,Corr,*Eav,*Cav,*E2av;
    bool ok;
    FILE *fp;
    char name[100],name_out[100];

    if(argc!=2) {printf("Usage: ./a.out [T]\n"); exit(0);}
    sscanf(argv[1],"%le", &T);
    init_random(0,0);
    N = 2048.;  // system size 
    conn = 5;   // connectivity
    NRRG = 8;   // number of realizations of the RRG
    ND = 8;     // number of dynamical trajectories
    kappa = 0.1;
    K1 = kappa;
    K2 = kappa/(double)conn;
    T_w = 16384;    // waiting time (initialization steps)
    quanti_out = (int) ((double)T_w/32.);   // number of output during Tw
    nmc = 14;   // sampling steps: N * 2^nmc
    quanti_punti = nmc + 5;
    Eav = (double *) calloc(sizeof(double), quanti_punti);
    E2av = (double *) calloc(sizeof(double), quanti_punti);
    Cav = (double *) calloc(sizeof(double), quanti_punti);
    time = (int *) calloc(sizeof(int), quanti_punti);
    time[0] = N/16;
    for(index=1;index<quanti_punti;index++)
    { 
        time[index] = time[index-1]*2;
        //printf("%d\n", time[index]);
    }
    for(i=0;i<quanti_punti;i++){ Eav[i] = 0.; E2av[i] = 0.; Cav[i] = 0.;}
    S = (int *) calloc(sizeof(int), N);         // spin configs
    S0 = (int *) calloc(sizeof(int), N);
    V = (int **) calloc(sizeof(int *), N);      // list of nearest neighbors
    quantiNN = (int *) calloc(sizeof(int), N);  // the number of NN neighbors for each site
    quantiNNN = (int *) calloc(sizeof(int), N); // the number of NNN neighbors for each site
    NN = (int **) calloc(sizeof(int *), N);     // list of NN neighbors
    NNN = (int **) calloc(sizeof(int *), N);    // list of NNN neighbors
    for(i=0;i<N;i++)
    {
        V[i] = (int *) calloc(sizeof(int), conn);
        NN[i] = (int *) calloc(sizeof(int), conn*(conn-1));
        NNN[i] = (int *) calloc(sizeof(int), conn*(conn-1)*(conn-1));
    }
    
    sprintf(name, "DYN-N%d-T%g-kappa%g-c%d", N, T, kappa, conn);
    sprintf(name_out, "out-N%d-T%g-kappa%g-c%d", N, T, kappa, conn);
    cnt = 0;
    for(r=0;r<NRRG;r++)
    {
        ok = false;
        while(!ok) ok = cRRG(conn,N,V);
        find_neighbors(conn,N,V,quantiNN,quantiNNN,NN,NNN);
        E = initialization(conn,N,V,quantiNN,quantiNNN,NN,NNN,K1,K2,S); // random distribution
        for(i=0;i<N;i++) S0[i] = S[i];
        Corr = 1.*N;
        for(i=0;i<T_w;i++)
        {
            for(j=0;j<N;j++) MCstep(conn,N,V,quantiNN,quantiNNN,NN,NNN,K1,K2,T,S,&E,S0,&Corr);
            if((i+1)%quanti_out==0)
            {
                fp = fopen(name_out, "a");
                fprintf(fp, "%d  %g  %g  %d\n", i, E, Corr, r);
                fclose(fp);
            }
        }
    
        for(j=0;j<ND;j++)
        {
            for(i=0;i<N;i++) S0[i] = S[i];
            Corr = 1.*N;
            i = 0;
            for(index=0;index<quanti_punti;index++){
                while(i<time[index]){
                    MCstep(conn,N,V,quantiNN,quantiNNN,NN,NNN,K1,K2,T,S,&E,S0,&Corr);
                    i++;
                }
                //printf("%d %d -> %d : %g %g\n", r, j, index, E/(double)N, Corr/(double)N);
                Eav[index] = (cnt*Eav[index] + E/(double)N)/(double)(cnt + 1);
                E2av[index] = (cnt*E2av[index] + E/(double)N*E/(double)N)/(double)(cnt + 1);
                Cav[index] = (cnt*Cav[index] + Corr/(double)N)/(double)(cnt + 1);
            }
            cnt++;
            fp = fopen(name, "w");
            for(index=0;index<quanti_punti;index++) fprintf(fp, "%g %g %g %g\n", time[index]/(double)N, Cav[index], Eav[index], sqrt(E2av[index] - Eav[index]*Eav[index]));
            fclose(fp);
        }
    }
    
    return 0;
}
