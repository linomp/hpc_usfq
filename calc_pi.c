#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define N 204800

double cuarto(double xi){
  return sqrt(1-pow(x,2));
}

int main(void) {

	float t_loop;

	struct timeval t_i, t_f;

	double dx = 1.0/(double)N;
	double cuartopi = 0;
	double val;

	gettimeofday(&t_i, NULL);

	for (int i = 0; i<N; i++){
    cuartopi += dx * cuarto(i*dx);
  }

	gettimeofday(&t_f, NULL);

	t_loop =  (double) t_f.tv_sec * 1000.0 + (double) t_f.tv_usec / 1000.0 - 
			 ((double) t_i.tv_sec * 1000.0 + (double) t_i.tv_usec / 1000.0);

	printf("Pi = %f \nTiempo del loop %f ms\n", 4 * cuartopi, t_loop );

	return 0;
}
