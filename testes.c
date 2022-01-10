#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void parameters(float *t0, int *N, float *rho, float *tmax, float *dt){
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
    fscanf(arq, "%f", rho);
    fscanf(arq, "%f", tmax);
    fscanf(arq, "%f", dt);

    fclose(arq);
}

void init(int N, float rho, float t0, float dt, 
            float *x, float *y, float *z,
            float *xm, float *ym, float *zm, 
            float *v_x, float *v_y, float *v_z){

    // Soma das velocidades
    float sumv_x = 0;
    float sumv_y = 0;
    float sumv_z = 0;

    // Soma das velocidades ao quadrado
    float sumv2_x = 0;
    float sumv2_y = 0;
    float sumv2_z = 0; 
    
    // Tamanho da caixa e espaçamento entre as partículas
    float L = cbrt(N/rho);
    int n3 = ceil(cbrt(N));
    float espaco = (float) L/(n3-1);

    printf("\nL: %f \nn3: %i \nespaco: %f\n", L, n3, espaco);

    int i = 0, j = 0, k = 0;
    for(int part = 0 ; part < N ; part++){
            
        x[part] = i * espaco;
        y[part] = j * espaco;
        z[part] = k * espaco;
    
        i++;
        if(i == n3){
            i = 0;
            j++;
            if(j == n3){
                j = 0;
                k++;
            }
        }
        
        v_x[part] = ((float)rand()/(float)(RAND_MAX) * 1) - 0,5;
        v_y[part] = ((float)rand()/(float)(RAND_MAX) * 1) - 0,5;
        v_z[part] = ((float)rand()/(float)(RAND_MAX) * 1) - 0,5;
    

        sumv_x = sumv_x + v_x[part];
        sumv_y = sumv_y + v_y[part];
        sumv_z = sumv_z + v_z[part];
        
        sumv2_x = sumv2_x + pow(v_x[part] , 2);
        sumv2_y = sumv2_y + pow(v_y[part] , 2);
        sumv2_z = sumv2_z + pow(v_z[part] , 2);    
    }
        

    // Media das velocidades
    sumv_x = sumv_x / N;
    sumv_y = sumv_y / N;
    sumv_z = sumv_z / N;
    
    // Media das velocidades quadradas
    sumv2_x = sumv2_x / N;
    sumv2_y = sumv2_y / N;
    sumv2_z = sumv2_z / N;

    // Fator de escala das velocidades
    float fs_x = sqrt(3 * t0 / sumv2_x);
    float fs_y = sqrt(3 * t0 / sumv2_y);
    float fs_z = sqrt(3 * t0 / sumv2_z);

    for(int i = 0 ; i < N ; i++){
        v_x[i] = (v_x[i] - sumv_x) * fs_x;
        v_y[i] = (v_y[i] - sumv_y) * fs_y;
        v_z[i] = (v_z[i] - sumv_z) * fs_z;
        
        
        xm[i] = x[i] - v_x[i] * dt;
        ym[i] = y[i] - v_y[i] * dt;
        zm[i] = z[i] - v_z[i] * dt;
    }
}

/*void force(int N, float *x, float box)
{
    float en = 0;
    float f[N] = {0};
    float xr;
    float r2;
    float rc2 = box/2;
    float r2i, r6i, ff; // variáveis para o cálculo da força
    float ecut = 4 * ((1/pow(rc2, 6)) - (1/pow(rc2, 3)) );
    
    for(int i = 0 ; i < N-1 ; i++)
    {
        for(int j = i+1 ; j < N ; j++)
        {
            xr = x[i] - x[j];
            xr = xr - box * round(xr/box);
            r2 = pow(xr, 2);

            if(r2 < rc2)
            {
                r2i = 1/r2;
                r6i = pow(r2i, 2);
                ff = 48 * r2i * r6i * (r6i - 0,5);
                f[i] = f[i] + ff * xr;
                f[j] = f[j] - ff * xr;

                en = en + 4 * r6i * (r6i - 1) - ecut;
            }
        }
    }



}*/

int main(){

    int N;
    float t0, rho, tmax, dt;
    
    parameters(&t0, &N, &rho, &tmax, &dt);
    printf(" t0: %f\n N: %d\n Density: %f\n tmax: %f\n dt: %f\n", t0, N, rho, tmax, dt);

    float x[N], y[N], z[N], xm[N], ym[N], zm[N], v_x[N], v_y[N], v_z[N];

    init(N, rho, t0, dt, x, y, z, xm, ym, zm, v_x, v_y, v_z);
    
    for(int i = 0 ; i<N; i++){
        printf("\nparticula %d\n", i);
        printf(" x: %f y: %f z: %f\n", x[i], y[i], z[i]);
        printf(" vx: %f vy: %f vz: %f\n", v_x[i], v_y[i], v_z[i]);
    }
}

