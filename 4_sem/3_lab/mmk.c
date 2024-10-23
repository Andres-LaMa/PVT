#define _POSIX_C_SOURCE 1
// #define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

# define M_PI 3.14159265358979323846

double getrand(unsigned int *seed){
    return (double)rand_r(seed)/RAND_MAX;
}

double function(double x, double y){
    return exp(2*x+2*y);
}

double s_mmk(const int n, double(*thisfunction)(double, double)){
    printf("Serial numerical integration by Monte Carlo method: n = %d\n", n);
    double time = omp_get_wtime();
    int in = 0;
    double s = 0;
    {
        double s_lock = 0;
        int in_lock = 0;
        unsigned int seed = omp_get_thread_num();
        for (int i = 0; i < n; i++){
            double x = getrand(&seed) * M_PI;
            double y = getrand(&seed);
            if (y <= 1-x){
                in_lock++;
                s_lock += thisfunction(x, y);
            }
        }
        s += s_lock;
        in += in_lock;
    }
    time = omp_get_wtime() - time;
    double v = M_PI * in / n;
    double res = v * s / in;
    printf ("%d %.12f\n", 1, time);
    return time;
}

void p_mmk(const int n, double(*thisfunction)(double, double), double serial_time){
    printf("Parallel numerical integration by Monte Carlo method: n = %d\n", n);
    for (int p = 2; p <= 8; p+=2)
    {
        double time = omp_get_wtime();
        int in = 0;
        double s = 0;
        #pragma omp parallel num_threads(p)
        {
            double s_lock = 0;
            int in_lock = 0;
            unsigned int seed = omp_get_thread_num();
            #pragma omp for nowait
            for (int i = 0; i < n; i++){
                double x = getrand(&seed) * M_PI;
                double y = getrand(&seed);
                if (y <= 1-x){
                    in_lock++;
                    s_lock += thisfunction(x, y);
                }
            }
            #pragma omp atomic
            s += s_lock;
            #pragma omp atomic
            in += in_lock;
        }
        time = serial_time/(omp_get_wtime() - time);
        // time = (omp_get_wtime() - time);
        double v = M_PI * in / n;
        double res = v * s / in;
        printf ("%d %.12f\n", p, time);
    }
}

int main(){
    printf("\n");
    for (int n = 1e7; n <= 1e8 ; n*=10){
        p_mmk(n, *function, s_mmk(n, *function));
        printf("\n");
    }
    return 0;
}