#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void parameters(float *t0, int *N, float *density, float *tmax, float *dt);
void init(int N, float t0, float dt, float *x, float *xm, float *v_x);
void force(f, en);
void integrate(f, en);
void sample();

int main(){

    int N;
    float t0, density, tmax, dt;
    parameters(&t0, &N, &density, &tmax, &dt); //atribui valores para as condicoes iniciais conforme arquivo de parametros do usuario

    float x[N], xm[N], v_x[N]; //define vetores para as posicoes e velocidades

    // Função realiza a inicialização das posições das partículas em uma rede cristalina com velocidades aleatórias
    init(N, t0, dt, x, xm, v_x);

    for(int t = 0; t < tmax; t += dt){
        force(f, en);
        integrate(f, en);
        sample();
    }

    return 0;
}

void parameters(float *t0, int *N, float *density, float *tmax, float *dt){
    FILE *arq;
    char nome[30];

    do{
        printf("Nome do arquivo de parametros: ");
        fflush(stdin);
        gets(nome);
        arq = fopen(nome, "r");
        if(arq == NULL)
            printf("Nao foi possivel abrir esse arquivo.\n");
    }while(arq == NULL);

    fscanf(arq, "%f", t0);
    fscanf(arq, "%d", N);
    fscanf(arq, "%f", density);
    fscanf(arq, "%f", tmax);
    fscanf(arq, "%f", dt);

    fclose(arq);
}


void init(int N, float t0, float dt, float *x, float *xm, float *v_x){

    float sumv = 0; // Soma das velocidades
    float sumv2 = 0; // Soma das velocidades ao quadrado

    float espaco = (float) 1/(N-1); // Espaçamento entre as partículas

    for(int i = 0 ; i < N ; i++){

        x[i] = i * espaco;
        v_x[i] = ((float)rand()/(float)(RAND_MAX) * 1) - 0,5;

        sumv = sumv + v_x[i];
        sumv2 = sumv2 + pow(v_x[i] , 2);
    }

    sumv = sumv / N;
    sumv2 = sumv2 / N;

    float fs = sqrt(3 * t0 / sumv2); // Fator de escala das velocidades

    for(int i = 0 ; i < N ; i++){
        v_x[i] = (v_x[i] - sumv) * fs;
        xm[i] = x[i] - v_x[i] * dt;
    }
}

void force(f, en){

}

void integrate(f, en){

}

void sample(){

}
