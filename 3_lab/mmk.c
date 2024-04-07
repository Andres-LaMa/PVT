// #define _POSIX_C_SOURCE 1
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

double getrand(unsigned int *seed){
    return (double)rand_r(seed)/RAND_MAX;
}

double function(double x, double y){
    return exp(x+y)*exp(x+y);
}

int main(){
    const int n = 10000000;
    printf("Numerical integration by Monte Carlo method: n = %d\n", n);
    int in = 0;
    double s = 0;
    #pragma omp parallel
    {
        double s_lock = 0;
        int in_lock = 0;
        unsigned int seed = omp_get_thread_num();
        #pragma omp for nowait
        for (int i = 0; i < n; i++){
            double x = getrand(&seed) * M_PI;
            double y = getrand(&seed);
            if (y <= sin(x)){
                in_lock++;
                s_lock += function(x, y);
            }
        }
        #pragma omp atomic
        s += s_lock;
        #pragma omp atomic
        in += in_lock;
    }
    double v = M_PI * in / n;
    double res = v * s / in;
    printf ("Result: %.12f, n %d\n", res, n);
    return 0;
}