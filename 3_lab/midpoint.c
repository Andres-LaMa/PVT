#include <stdio.h>
#include <math.h>
#include <omp.h>

double function(double x){
    // return x/(sin(2*x)*sin(2*x)*sin(2*x));
    return (-x*x);
}

double s_midpoint_method(double a, double b, double(*thisfunction)(double)){
    const double eps = 1e-6;
    const int n0 = 1e9;
    double t = omp_get_wtime();
    printf("\nSerial numerical integration: [%f, %f], n0 = %d, EPS = %f\n", a, b, n0, eps);
    double sq[2];
    int n = n0, k;
    double delta = 1;
    for (k = 0; delta > eps; n*=2, k^=1){
        double h = (b-a)/n;
        double s = 0.0;
        for (int i = 0; i < n; i++)
            s+=thisfunction(a+h*(i+0.5));
        sq[k] += s*h;
        if (n>n0)
            delta = fabs(sq[k] - sq[k^1])/3.0;
        #if 0
        printf("n=%d i=%d sq=%.12f delta=%.12f\n", n, k, sq[k], delta);
        #endif
    }
    #if 0
    printf("Result Pi: %.12f, Runge rule: EPS %e, n %d\n", sq[k]*sq[k], eps, n/2);
    #endif
        
    t = (omp_get_wtime() - t);
    printf("%d %.6f\n", 1, t);
    return t;
}

void p_midpoint_method(double a, double b, double(*thisfunction)(double), double serial_time){
    const double eps = 1e-6;
    const int n0 = 1e9;
    printf("\nParallel numerical integration: [%f, %f], n0 = %d, EPS = %f\n", a, b, n0, eps);
    for (int p = 2; p <= 8; p+=2){
        double t = omp_get_wtime();
        double sq[2] = {0};
        #pragma omp parallel num_threads(p)
        {
            int n = n0, k;
            double delta = 1;
            for (k = 0; delta > eps; n*=2, k^=1){
                double h = (b-a)/n;
                double s = 0.0;
                #pragma omp barrier
                #pragma omp for nowait
                for (int i = 0; i < n; i++)
                    s+=thisfunction(a+h*(i+0.5));
                #pragma omp atomic
                sq[k] += s*h;
                #pragma omp barrier
                if (n>n0)
                    delta = fabs(sq[k] - sq[k^1])/3.0;
                #if 0
                printf("n=%d i=%d sq=%.12f delta=%.12f\n", n, k, sq[k], delta);
                #endif
            }
            // #pragma omp master
            // printf("Result: %.12f, Runge rule: EPS %e, n %d\n", sq[k]*sq[k], eps, n/2);
        }
        t = serial_time/(omp_get_wtime() - t);
        printf("%d %.6f\n", p, t);
    }
}

int main(){
    double a = 0.1, b = 0.5;
    p_midpoint_method(a, b, *function, s_midpoint_method(a, b, *function));
    
    return 0;
}