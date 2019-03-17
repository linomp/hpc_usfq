// 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h> 

#define N 10000000
void printSum(float* sum, int length){
  for (int i = 0; i < length; i++) {
		printf("%.2f, ", sum[i]); 
  }
  printf("\n");
}

void initializeVect(float* v, float val, int length){
	for (int i = 0; i < length; i++) {
		v[i] = val; 
  }
}

int main(int argc, char*argv[]) {
	float t_serial, t_par, t_total;
	struct timeval t_ii, t_i, t_s, t_f; 
	
	float *v1, *v2, *res;

	gettimeofday(&t_ii, NULL); 

	int nthreads = atoi(argv[1]);

	int length = N;
  int size = length * sizeof(float);
 
  v1 = (float *)malloc(size);
  v2 = (float *)malloc(size);
	result = (float *)malloc(size);

	initializeVect(v1, 50, length);
	initializeVect(v2, 30, length);

	gettimeofday(&t_i, NULL); 		 
	for (int i = 0; i<length; i++){ 
		result[i] = v1[i] + v2[i];
	}
	gettimeofday(&t_s, NULL);
	t_serial =  (double) t_s.tv_sec * 1000.0 + (double) t_s.tv_usec / 1000.0 - 
			 ((double) t_i.tv_sec * 1000.0 + (double) t_i.tv_usec / 1000.0);
	
	gettimeofday(&t_i, NULL); 	
	#pragma omp parallel for num_threads(nthreads)
	for (int i = 0; i<length; i++){ 
		result[i] = v1[i] + v2[i];
	}	 
	gettimeofday(&t_f, NULL);
  t_par =  (double) t_f.tv_sec * 1000.0 + (double) t_f.tv_usec / 1000.0 - 
			 ((double) t_i.tv_sec * 1000.0 + (double) t_i.tv_usec / 1000.0);
		
	t_total =  (double) t_f.tv_sec * 1000.0 + (double) t_f.tv_usec / 1000.0 - 
			 ((double) t_ii.tv_sec * 1000.0 + (double) t_ii.tv_usec / 1000.0);

	float speedup = t_serial/t_total;
	float eff = speedup/nthreads;

	//printSum(result, length); 
	printf("\nHilos: %i\nSpeedup: %f\nEficiencia: %f\nTiempo serial: %f ms\nTiempo paralelo: %f ms\nTiempo total %f ms\n", nthreads, speedup, eff, t_serial, t_par, t_total );

	free(v1);
	free(v2);
	free(result); 

	return 0;
}