#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>
#include <omp.h>

int N = 25000;
int M = 25000;

double wtime()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1E-9;
}

void matrix_serial(double *a, double *b, double *c, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        c[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            c[i] += a[i * n + j] * b[j];
        }
    }
}

void run_serial(double *time)
{
    double *a, *b, *c;
    a = malloc(sizeof(*a) * N * M);
    b = malloc(sizeof(*b) * N);
    c = malloc(sizeof(*c) * M);

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            a[i * N + j] = i + j;
        }
    }
    for (int i = 0; i < N; i++)
    {
        b[i] = i;
    }
    double t = wtime();
    matrix_serial(a, b, c, M, N);
    t = (*time) = wtime() - t;
    printf("Elapsed time (serial): %.6f sec.\n", t);
    free(a);
    free(b);
    free(c);
}

void matrix_parallel(double *a, double *b, double *c, int m, int n)
{
#pragma omp parallel // num_threads(2)
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);

        for (int i = lb; i <= ub; i++)
        {
            c[i] = 0.0;
            for (int j = 0; j < n; j++)
            {
                c[i] += a[i * n + j] * b[j];
            }
        }
    }
}

void run_parallel(double *s_time)
{
    double *a, *b, *c;
    a = malloc(sizeof(*a) * N * M);
    b = malloc(sizeof(*b) * N);
    c = malloc(sizeof(*c) * M);

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            a[i * N + j] = i + j;
        }
    }
    for (int j = 0; j < N; j++)
        b[j] = j;

    double t = wtime();
    matrix_parallel(a, b, c, M, N);
    t = (*s_time) = wtime() - t;
    printf("Elapsed time (parallel): %.6f sec.\n", t);
    free(a);
    free(b);
    free(c);
}

int main()
{
    FILE *f = fopen("data", "w");
    double p_time, s_time;
    for (int j = 15000; j <= 25000; j += 5000)
    {
        for (int i = 0; i <= 8; i += 2)
        {
            M = N = j;
            omp_set_num_threads(i);
            printf("Memory used: %" PRIu64 " Mib\n", ((M * N + M + N) * sizeof(double) >> 20));
            run_serial(&s_time);
            run_parallel(&p_time);
            printf("Boost algoritmic = %lf on %d threads\n\n", s_time / p_time, i);
            fprintf(f, "%d\t%d\t%lf\n", j, (i != 0 ? i : 1), s_time / p_time);
        }
    }
    fclose(f);
}