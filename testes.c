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

void init(int N, float rho, float t0, float dt, float *L,
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
    *L = cbrt(N/rho);
    int n3 = ceil(cbrt(N));
    float espaco = (float) *L/(n3);

    printf("\nL: %f \nn3: %i \nespaco: %f\n", *L, n3, espaco);

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

void force(int N, double *fx, double *fy, double *fz, double *en, float *x, float *y, float *z, float L)
{   
    *en = 0;

    for(int i = 0 ; i < N ; i++){
        fx[i] = fy[i] = fz[i] = 0;
    }

    // Declaração das variáveis para o cálculo da força
    float xr, yr, zr;
    float r2, r2i, r6i, ff; 
    
    // Distância máxima de interação e potencial de Lennard Jones
    float rc2 = pow(L/2, 2); 
    printf("rc2 = %f\n", rc2);
    float ecut = 4 * ((1/pow(rc2, 6)) - (1/pow(rc2, 3)));
    
    for(int i = 0 ; i < N-1 ; i++)
    {
        for(int j = i+1 ; j < N ; j++)
        {
            // Cálculo das distâncias de um par de partículas em cada coordenada
            xr = x[i] - x[j];
            xr = xr - L * round(xr/L);

            yr = y[i] - y[j];
            yr = yr - L * round(yr/L);
            
            zr = z[i] - z[j];
            zr = zr - L * round(zr/L);
            
            // Cálculo da distância total de um par de partículas
            r2 = pow(xr, 2) + pow(yr, 2) + pow(zr, 2);
            printf("r2 particula %d = %f\n", i, r2);

            if(r2 <= rc2)
            {
                r2i = 1/r2;
                r6i = pow(r2i, 3);
                ff = 48 * r2i * r6i * (r6i - 0,5);
                printf("ff particula %d = %f\n", i, ff);
                
                fx[i] = fx[i] - ff * xr;
                fx[j] = fx[j] - ff * xr;
                
                fy[i] = fy[i] + ff * yr;
                fy[j] = fy[j] - ff * yr;
                
                fz[i] = fz[i] + ff * zr;
                fz[j] = fz[j] - ff * zr;
    
                *en = *en + 4 * r6i * (r6i - 1) - ecut;
            }
        }
    }
}

void integrate(float N, float dt, double en,
                double *fx, double *fy, double *fz,
                float *x, float *y, float *z,
                float *xm, float *ym, float *zm,
                float *temp, float *etot){

    float sumv = 0, sumv2 = 0, xx, yy, zz, vi, vj, vk;

    for(int i = 0 ; i < N ; i++){

        xx = 2*x[i] - xm[i] + pow(dt, 2) * fx[i];
        vi = (xx - xm[i]) / (2*dt);

        yy = 2*y[i] - ym[i] + pow(dt, 2) * fy[i];
        vj = (yy - ym[i]) / (2*dt);

        zz = 2*z[i] - zm[i] + pow(dt, 2) * fz[i];
        vk = (zz - zm[i]) / (2*dt);

        sumv = sumv + vi + vj + vk;
        sumv2 = sumv2 + pow(vi, 2) + pow(vj, 2) + pow(vk, 2);
        
        xm[i] = x[i];
        x[i] = xx;

        ym[i] = y[i];
        y[i] = yy;

        zm[i] = z[i];
        z[i] = zz;
    }
    
    if(sumv != 0){
        printf("ERRO: A velocidade do centro de massa é diferente de zero\n");
        //Break
    }

    *temp = sumv2/(3*N);
    *etot = (en + 0.5*sumv2)/N;
}


int main(){

    int N;
    float t0, rho, tmax, dt, L;
    float x[N], y[N], z[N], xm[N], ym[N], zm[N], v_x[N], v_y[N], v_z[N];
    double en, fx[N], fy[N], fz[N];
    float temp, etot;

    parameters(&t0, &N, &rho, &tmax, &dt);
    printf(" t0: %f\n N: %d\n Density: %f\n tmax: %f\n dt: %f\n", t0, N, rho, tmax, dt);

    init(N, rho, t0, dt, &L, x, y, z, xm, ym, zm, v_x, v_y, v_z);
    force(N, fx, fy, fz, &en, x, y, z, L);
    integrate(N, dt, en, fx, fy, fz, x, y, z, xm, ym, zm, &temp, &etot);
    
    for(int i = 0 ; i<N; i++){
        printf("\nParticula %d\n", i);
        printf(" x: %f y: %f z: %f\n", x[i], y[i], z[i]);
        printf(" vx: %f vy: %f vz: %f\n", v_x[i], v_y[i], v_z[i]);
        printf(" fx: %f fy: %f fz: %f\n", fx[i], fy[i], fz[i]);
    }
    
    printf("Energia Total: %f", en);
    printf("Energia por Partícula: %f", etot);
    printf("Temperatura: %f", temp);

    return 0;
}
