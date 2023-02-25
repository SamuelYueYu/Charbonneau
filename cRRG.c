#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>

extern void init_random(int seed,int seed2);
extern double Xrandom(void);

bool cRRG(int conn, int N, int **V){
    int cnt,i,j,v1,v2,s,r,u1,u2,*U,*Nv;
    int *siti_restanti,cnt_siti,**Link,cnt_link,l,cnt_link_aux,**Link_aux,nuovi_link;
    bool ok,estrai,riuscito;
    double p,prob;

    U = (int *) calloc(sizeof(int), conn*N);
    Nv = (int *) calloc(sizeof(int), N);
    siti_restanti = (int *) calloc(sizeof(int), 2*conn);
    Link =  (int **) calloc(sizeof(int *), conn*(2*conn - 1));
    Link_aux =  (int **) calloc(sizeof(int *), conn*(2*conn - 1));
    for(i=0;i<conn*(2*conn - 1);i++){
        Link[i] = (int *) calloc(sizeof(int), 2);
        Link_aux[i] = (int *) calloc(sizeof(int), 2);
    }
    
    riuscito = true;
    for(i=0;i<N;i++){
        Nv[i] = 0;
        for(j=0;j<conn;j++){
            U[conn*i+j] = i;
        }
    }
    cnt = conn*N;
    while(cnt>2*conn){
        estrai = false;
        while(!estrai){
            estrai = true;
            v1 = (int) (cnt*Xrandom());
            s = U[v1];
            v2 = (int) (cnt*Xrandom());
            r = U[v2];
            if(s==r) estrai = false;
            else{
                for(i=0;i<Nv[s]&&estrai;i++){
                    if(V[s][i]==r) estrai = false;
                }
            }
        }

        V[s][Nv[s]] = r;
        V[r][Nv[r]] = s;
        Nv[s]++;
        Nv[r]++;
        if(v1<v2){
            u1 = v1;
            u2 = v2;
        }
        else{
            u1 = v2;
            u2 = v1;
        }
        if(u2==cnt-2) U[u1] = U[cnt-1];
        else{
            if(u2==cnt-1){
                if(u1<cnt-2) U[u1] = U[cnt-2];
            }
            else{
                U[u2] = U[cnt-1];
                U[u1] = U[cnt-2];
            }
        }
        cnt = cnt - 2;
    }

    //  FASE 2
    siti_restanti[0] = U[0];
    cnt_siti = 1;
    for(i=1;i<2*conn;i++){
        ok = true;
        for(j=0;j<cnt_siti&&ok;j++) if(U[i]==siti_restanti[j]) ok = false;
        if(ok){
            siti_restanti[cnt_siti] = U[i];
            cnt_siti++;
        }
    }
    cnt_link = 0;
    for(i=0;i<cnt_siti;i++){
        for(j=i+1;j<cnt_siti;j++){
            v1 = siti_restanti[i];
            v2 = siti_restanti[j];
            ok = true;
            for(l=0;l<Nv[v1]&&ok;l++){
                if(V[v1][l]==v2) ok = false;
            }
            if(ok){
                Link[cnt_link][0] = v1;
                Link[cnt_link][1] = v2;
                cnt_link++;
            }
        }
    }
    nuovi_link = 0;
    while(nuovi_link<conn&&riuscito){
        l = (int) (cnt_link*Xrandom());
        v1 = Link[l][0];
        v2 = Link[l][1];
        p = Xrandom();
        prob = ((double)(conn-Nv[v1]))*((double)(conn-Nv[v2]))/(double)(conn*conn);
        if(p<prob){
            V[v1][Nv[v1]] = v2;
            V[v2][Nv[v2]] = v1;
            Nv[v1]++;
            Nv[v2]++;
            nuovi_link++;
            cnt_link--;
            for(j=l;j<cnt_link;j++){
                Link[j][0] = Link[j+1][0];
                Link[j][1] = Link[j+1][1];
            }
            cnt_link_aux = 0;
            for(j=0;j<cnt_link;j++){
                v1 = Link[j][0];
                v2 = Link[j][1];
                if(Nv[v1]<conn&&Nv[v2]<conn){
                    Link_aux[cnt_link_aux][0] = v1;
                    Link_aux[cnt_link_aux][1] = v2;
                    cnt_link_aux++;
                }
            }
            for(j=0;j<cnt_link_aux;j++){
                Link[j][0] = Link_aux[j][0];
                Link[j][1] = Link_aux[j][1];
            }
            cnt_link = cnt_link_aux;
            if(cnt_link==0&&nuovi_link<conn) riuscito = false;
        }
    }

    free(U);
    free(Nv);
    for(i=0;i<conn*(2*conn - 1);i++){
        free(Link[i]);
        free(Link_aux[i]);
    }
    free(Link);
    free(Link_aux);
    
    return riuscito;
}

void find_neighbors(int conn, int N, int **V, int *quantiNN, int *quantiNNN, int **NN, int **NNN){
    int i,s1,s2,sa,v,v2;
    bool diverso1,diverso2,diverso3;
    
    for(i=0;i<N;i++){
        quantiNN[i] = 0;
        for(s1=0;s1<conn;s1++){
            v = V[i][s1];
            for(s2=0;s2<conn;s2++){
                v2 = V[v][s2];
                if(v2!=i){
                    diverso1 = true;
                    for(sa=0;sa<conn&&diverso1;sa++){
                        if(V[i][sa]==v2) diverso1 = false;
                    }
                    if(diverso1){
                        diverso2 = true;
                        for(sa=0;sa<quantiNN[i]&&diverso2;sa++){
                            if(NN[i][sa]==v2) diverso2 = false;
                        }
                        if(diverso2){
                            NN[i][quantiNN[i]] = v2;
                            quantiNN[i]++;
                        }
                    }
                }
            }
        }
    }
    for(i=0;i<N;i++){
        quantiNNN[i] = 0;
        for(s1=0;s1<quantiNN[i];s1++){
            v = NN[i][s1];
            for(s2=0;s2<conn;s2++){
                v2 = V[v][s2];
                if(v2!=i){
                    diverso1 = true;
                    for(sa=0;sa<conn&&diverso1;sa++){
                        if(V[i][sa]==v2) diverso1 = false;
                    }
                    if(diverso1){
                        diverso2 = true;
                        for(sa=0;sa<quantiNN[i]&&diverso2;sa++){
                            if(NN[i][sa]==v2) diverso2 = false;
                        }
                        if(diverso2){
                            diverso3 = true;
                            for(sa=0;sa<quantiNNN[i]&&diverso3;sa++){
                                if(NNN[i][sa]==v2) diverso3 = false;
                            }
                            if(diverso3){
                                NNN[i][quantiNNN[i]] = v2;
                                quantiNNN[i]++;
                            }
                        }
                    }
                }
            }
        }
    }
}
