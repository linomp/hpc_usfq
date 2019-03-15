#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define N 204800

double fcn_cuarto(double xi){
	return sqrt(1-pow(xi, 2));
}

int main(void) {

	float t_loop;

	struct timeval t_i, t_f;

	double dx = 1.0/(double)N;
	double cuartopi = 0;
	double val;

	gettimeofday(&t_i, NULL);

	// hace y corre 2 copias (INCORRECTO)
	/*#pragma omp parallel
	{
		for (int i = 0; i<N-1; i++){ 
			val = dx*fcn_cuarto(i*dx);
			#pragma omp critical
				cuartopi += val;
		}
	}*/

	/*	
	#pragma omp parallel for
	for (int i = 0; i<N-1; i++){ 
		val = dx*fcn_cuarto(i*dx);
		#pragma omp critical
			cuartopi += val;
	}
	*/

		
	#pragma omp parallel for reduction(+: cuartopi)
	for (int i = 0; i<N-1; i++){ 
		cuartopi += dx*fcn_cuarto(i*dx);
	}
	
	

	gettimeofday(&t_f, NULL);

	t_loop =  (double) t_f.tv_sec * 1000.0 + (double) t_f.tv_usec / 1000.0 - 
			 ((double) t_i.tv_sec * 1000.0 + (double) t_i.tv_usec / 1000.0);

	printf("Pi = %f \nTiempo del loop %f ms\n", 4 * cuartopi, t_loop );

	return 0;
}