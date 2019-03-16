#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

//#define N 2048000

int main(int argc, char** argv)
{
  
  long N = 2048000*(int) argv[0];

  float t_loop, t_total;
  struct timeval t_i, t_f, t_fh;

  double dx = 1.0 / (double)N;
  double cuartopi = 0;
  int iters = 1000000;

  gettimeofday(&t_i, NULL);
  #pragma omp parallel for reduction(+ \
                                   : cuartopi)
  for (int i = 0; i < N; ++i)
  {
    cuartopi += dx * sqrt(1 - pow(i * dx, 2));
  }
  gettimeofday(&t_f, NULL);

  double h = 2.234;

  for (int i = 0; i < iters; ++i)
  {
    h = h * cos(cuartopi * h);
  }
  gettimeofday(&t_fh, NULL);

  t_loop = (double)t_f.tv_sec * 1000.0 + (double)t_f.tv_usec / 1000.0 -
           ((double)t_i.tv_sec * 1000.0 + (double)t_i.tv_usec / 1000.0);

  t_total = (double)t_fh.tv_sec * 1000.0 + (double)t_fh.tv_usec / 1000.0 -
            ((double)t_i.tv_sec * 1000.0 + (double)t_i.tv_usec / 1000.0);

  printf("Pi = %f\n h = %f\nTiempo parte paralela %f ms\nTiempo parte serial %f ms\nTiempo total % f ms\n ", (4 * cuartopi), h, t_loop,t_total-t_loop, t_total);

  return 0;
}
