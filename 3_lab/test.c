#include <stdio.h>
#include <math.h>
#include <omp.h>
double func(double x)
{
  return (sqrt(x * (3 - x)) / (x + 1));
}

void runge_serial(double eps, double a, double b, int n0)
{
  double sq[2];
  int n = n0, k;
  double delta = 1;
  for (k = 0; delta > eps; n *= 2, k ^= 1)
  {
    double h = (b - a) / n;
    double s = 0.0;
    sq[k] = 0;
    for (int i = 0; i < n; i++)
      s += func(a + h * (i + 0.5));
    sq[k] += s * h;
    if (n > n0)
      delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
    #if 0
    printf("n=%d i=%d sq=%.12f delta=%.12f\n", n, k, sq[k], delta);
    #endif
  }
  printf("Result of serial: %.12f; Runge rule: EPS %e, n %d\n", sq[k] * sq[k], eps, n / 2);
}

void runge_parallel(double eps, double a, double b, int n0, int count_threads)
{
  double sq[2];
  #pragma omp parallel num_threads(count_threads)
  {
    int n = n0, k;
    double delta = 1;
    for (k = 0; delta > eps; n *= 2, k ^= 1) {
      double h = (b - a) / n;
      double s = 0.0;
      sq[k] = 0;
      // Ждем пока все потоки закончат обнуление sq[k], s
      #pragma omp barrier
      #pragma omp for nowait
      for (int i = 0; i < n; i++)
        s += func(a + h * (i + 0.5));
      #pragma omp atomic
      sq[k] += s * h;
      // Ждем пока все потоки обновят sq[k]
      #pragma omp barrier
      if (n > n0)
        delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
      #if 0
      printf("n=%d i=%d sq=%.12f delta=%.12f\n", n, k, sq[k], delta);
    #endif
    }
    #pragma omp master
    printf("Result of %d threads: %.12f; Runge rule: EPS %e, n %d\n", count_threads, sq[k] * sq[k], eps, n / 2);
  }
}

int main(int argc, char **argv)
{
  double t_parallel = 0, t_serial = 0;
  const double eps = 1E-5;
  const double a = 1;
  const double b = 1.2;
  const int n0 = 100000000;
  double boost = 0;
  
  printf("Numerical integration: [%f, %f], n0 = %d, EPS = %f\n", a, b, n0, eps);
  t_serial = omp_get_wtime();
  runge_serial(eps, a, b, n0);
  t_serial = omp_get_wtime() - t_serial;
  
  printf("Serial time: %.6f\n", t_serial);
  for(int i = 2; i <= 8; i += 2)
  {
    t_parallel = omp_get_wtime();
    runge_parallel(eps, a, b, n0, i);
    t_parallel = omp_get_wtime() - t_parallel;
    boost = t_serial / t_parallel;
    
    printf("Parallel time of %d threads: %.6f\n", i, t_parallel);
  }
  // printf("Elapsed time (sec.): %.6f\n", t);
  return 0;
}