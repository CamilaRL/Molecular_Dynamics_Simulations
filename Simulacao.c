#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define t0 0.728
#define N 108
#define rho 0.8442
#define dt 0.001
#define tmax 1/dt


void init(float *L, float *x, float *y, float *z,
            float *xm, float *ym, float *zm,
            float *v_x, float *v_y, float *v_z);
void force(float *fx, float *fy, float *fz,
            double *U, float *x, float *y, float *z, float L);
void integrate(double U, float *temp, double *etot, double *k, float L,
                float *fx, float *fy, float *fz,
                float *x, float *y, float *z,
                float *xm, float *ym, float *zm);

int main(){

    // Declaracao de variaveis
    float x[N], y[N], z[N];
    float xm[N], ym[N], zm[N];
    float v_x[N], v_y[N], v_z[N];
    float fx[N], fy[N], fz[N];
    float L, temp;
    double U = 0, etot = 0, k = 0;

    FILE *arq;
    arq = fopen("saida.txt", "w");
    fprintf(arq, "# U K E\n");

    printf("\n\tCondicoes Iniciais\n");
    printf("\tt0: %f\n\tN: %d\n\trho: %d\n\ttmax: %d\n\tdt: %d\n", t0, N, rho, tmax, dt);

    init(&L, x, y, z, xm, ym, zm, v_x, v_y, v_z);
    
    for(int i = 0 ; i < tmax ; i++)
    {
        force(fx, fy, fz, &U, x, y, z, L);
        integrate(U, &temp, &etot, &k , L, fx, fy, fz, x, y, z, xm, ym, zm);
        fprintf(arq, "%f %f %f\n", U, k, etot);
    }
    
    fclose(arq);    
/*
    for(int i = 0 ; i<N; i++){
        printf("\nParticula %d\n", i);
        printf(" x: %f y: %f z: %f\n", x[i], y[i], z[i]);
        printf(" vx: %f vy: %f vz: %f\n", v_x[i], v_y[i], v_z[i]);
        printf(" fx: %f fy: %f fz: %f\n", fx[i], fy[i], fz[i]);
    }
*/
    printf("\nEnergia Potencial Total: %f", U);
    printf("\nEnergia por Particula: %f", etot);
    printf("\nTemperatura: %f\n", temp);

    return 0;
}

void init(float *L, float *x, float *y, float *z,
            float *xm, float *ym, float *zm,
            float *v_x, float *v_y, float *v_z){
    
    /*
        Funcao para posicionar as particulas em uma rede cristalina e
        designar a velocidade inicial de cada uma aleatoriamente.

        Entrada:
            L = tamanho da lateral da caixa
            x, y, z = arrays com as posicoes de cada particula em cada coordenada
            xm, ym, zm = arrays com as posicoes anteriores em cada coordenada
            v_x, v_y, vz = arrays com as velocidades em cada coordenada
    */

    // Definicao de variáveis
    float espaco, fs_x, fs_y, fs_z;
    int n3, part, i = 0, j = 0, k = 0;
    float sumv_x = 0, sumv_y = 0, sumv_z = 0;
    float sumv2_x = 0, sumv2_y = 0, sumv2_z = 0;

    // Tamanho da caixa e espacamento entre as particulas
    *L = cbrt((float)N/rho);
    n3 = ceil(cbrt(N));
    espaco = (float) *L/n3;

    printf("\n\tCaracteristicas da caixa\n");
    printf("\tL: %f\n\tn3: %i\n\tespaco: %f\n", *L, n3, espaco);

    // Inicializacao de cada particula com posicoes e velocidades
    for(part = 0 ; part < N ; part++){

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

        // Velocidades aleatorias
        v_x[part] = ((float)rand()/(float)(RAND_MAX) * 1) - 0.5;
        v_y[part] = ((float)rand()/(float)(RAND_MAX) * 1) - 0.5;
        v_z[part] = ((float)rand()/(float)(RAND_MAX) * 1) - 0.5;

        sumv_x += v_x[part];
        sumv_y += v_y[part];
        sumv_z += v_z[part];

        sumv2_x += pow(v_x[part], 2);
        sumv2_y += pow(v_y[part], 2);
        sumv2_z += pow(v_z[part], 2);
    }

    // Media das velocidades
    sumv_x /= N;
    sumv_y /= N;
    sumv_z /= N;

    // Media das velocidades quadradas
    sumv2_x /= N;
    sumv2_y /= N;
    sumv2_z /= N;

    // Fator de escala das velocidades
    fs_x = sqrt(3 * t0 / sumv2_x);
    fs_y = sqrt(3 * t0 / sumv2_y);
    fs_z = sqrt(3 * t0 / sumv2_z);

    // Correcao das velocidades
    for(i = 0 ; i < N ; i++){

        v_x[i] = (v_x[i] - sumv_x) * fs_x;
        v_y[i] = (v_y[i] - sumv_y) * fs_y;
        v_z[i] = (v_z[i] - sumv_z) * fs_z;

        // Posicao anterior conforme a conservacao de momentum
        xm[i] = x[i] - (v_x[i] * dt);
        ym[i] = y[i] - (v_y[i] * dt);
        zm[i] = z[i] - (v_z[i] * dt);
    }

}

void force(float *fx, float *fy, float *fz,
            double *U, float *x, float *y, float *z, float L){
    
    /*
        Funcao para calcular as forcas entre cada particula a partir do potencial de Lennard-Jones 
        e a energia potencial.

        Entrada:
            L = tamanho da lateral da caixa
            x, y, z = arrays com as posicoes de cada particula em cada coordenada
            fx, fy, fz = arrays com as forcas em cada coordenada
            en = energia potencial total
    */

    // Definição de variaveis
    float xr, yr, zr;
    float rc2, ecut;
    float r2, r2i, r6i, ff;

    // Forcas e energias nulas
    double en = 0;

    for(int i = 0 ; i < N ; i++){
        fx[i] = fy[i] = fz[i] = 0;
    }

    // Distância máxima de interação e potencial de Lennard-Jones
    rc2 = pow(L/2, 2);
    ecut = 4 * ((1/pow(rc2, 6)) - (1/pow(rc2, 3)));

    // Calculo das forcas para cada par de particulas
    for(int i = 0 ; i < N-1 ; i++)
    {
        for(int j = i+1 ; j < N ; j++)
        {
            // Cálculo das distâncias de um par de partículas em cada coordenada
            xr = x[i] - x[j];
            xr = xr - (L * round(xr/L));

            yr = y[i] - y[j];
            yr = yr - (L * round(yr/L));

            zr = z[i] - z[j];
            zr = zr - (L * round(zr/L));

            // Cálculo da distância total de um par de partículas
            r2 = pow(xr, 2) + pow(yr, 2) + pow(zr, 2);
            
            // Calculo da forca se r2 for menor que a distancia mínima
            if(r2 < rc2 && r2 != 0)
            {
                r2i = 1/r2;
                r6i = pow(r2i, 3);
                ff = 48 * r2i * r6i * (r6i - 0.5);

                fx[i] = fx[i] + ff * xr;
                fx[j] = fx[j] - ff * xr;

                fy[i] = fy[i] + ff * yr;
                fy[j] = fy[j] - ff * yr;

                fz[i] = fz[i] + ff * zr;
                fz[j] = fz[j] - ff * zr;

                // Atualizacao da energia potencial total
                en += (4 * r6i * (r6i - 1)) - ecut;
            }
            r2 = 0.0;
        }
    }
    *U = en/N;
}

void integrate(double U, float *temp, double *etot, double *k, float L,
                float *fx, float *fy, float *fz,
                float *x, float *y, float *z,
                float *xm, float *ym, float *zm){

    /*
        Funcao que computa as proximas posicoes e velocidades.

        Entrada:
            en = energia potencial total
            temp = temperatura do sistema
            etot = energia total do sistema (potencial + cinetica)
            fx, fy, fz = arrays com as forcas em cada coordenada
            x, y, z = arrays com as posicoes de cada particula em cada coordenada
            xm, ym, zm = arrays com as posicoes anteriores em cada coordenada
    */

    // Definicao de variaveis
    float sumvi = 0, sumvj = 0, sumvk = 0, sumv2 = 0;
    float xx = 0, yy = 0, zz = 0;
    float vi, vj, vk;

    for(int i = 0 ; i < N ; i++){

        // Calculo das proximas posicoes e velocidades
        xx = (2*x[i]) - xm[i] + (pow(dt, 2) * fx[i]);
        //printf("%f", xx);
        xx = xx - (L * round(xx/L));
        //printf("    %f\n", xx);
        vi = (xx - xm[i]) / (2*dt);
        
        yy = (2*y[i]) - ym[i] + (pow(dt, 2) * fy[i]);
        yy = yy - (L * round(yy/L));
        vj = (yy - ym[i]) / (2*dt);

        zz = (2*z[i]) - zm[i] + (pow(dt, 2) * fz[i]);
        zz = zz - (L * round(zz/L));
        vk = (zz - zm[i]) / (2*dt);

        sumvi += vi;
        sumvj += vj;
        sumvk += vk;
        sumv2 += pow(vi, 2) + pow(vj, 2) + pow(vk, 2);

        // Atualizacao das posicoes
        xm[i] = x[i];
        x[i] = xx;

        ym[i] = y[i];
        y[i] = yy;

        zm[i] = z[i];
        z[i] = zz;
    }

    // O centro de massa nao pode se mover
    if(-0.0001 < sumvi < 0.0001 || -0.0001 < sumvj < 0.0001 || -0.0001 < sumvk < 0.0001)
       printf("\nATENCAO: A velocidade do centro de massa e diferente de zero.\nVi = %f Vj = %f Vk = %f\n", sumvi, sumvj, sumvk);

    // Calculo da temperatura e energia total do sistema
    *temp = sumv2/(3*N);
    *k = 0.5*sumv2/N;
    *etot = U + *k;
}