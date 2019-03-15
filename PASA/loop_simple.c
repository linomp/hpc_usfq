#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define BENCH_REPEAT 10000

#define N 20480

void random_ints(int *a);
void random_floats(float *a);
int checkResult(int *a, int *b, int *c, int n);

int main(void) {

	float *a, *b, *c;
	unsigned int size = N * sizeof(int);

	float coef = 2.234;

	float t_loop;

	a = (float *)malloc(size); random_floats(a);
	b = (float *)malloc(size); random_floats(b);
	c = (float *)malloc(size);


	struct timeval t_i, t_f;

	gettimeofday(&t_i, NULL);

	#pragma omp parallel num_threads(4)
	{
	for (int r = 0; r < BENCH_REPEAT; ++r)
	{
	
		for (int i = 0; i < N; ++i)
		{
			c[i] += a[i] / coef + b[i];
		} 

	}
	}
	gettimeofday(&t_f, NULL);

	t_loop = (double) t_f.tv_sec * 1000.0 + (double) t_f.tv_usec / 1000.0 - 
			 ((double) t_i.tv_sec * 1000.0 + (double) t_i.tv_usec / 1000.0);

	printf("Tiempo del loop %f ms\n", t_loop );

	return 0;
}

void random_ints(int *a)
{
	for (unsigned int i = 0; i < N; i++){
		a[i] = rand();
	}
}

void random_floats(float *a)
{
	for (unsigned int i = 0; i < N; i++){
		a[i] = 1.0f / (float) rand();
	}
}

int checkResult(int *a, int *b, int *c, int n){
	int correcto = 1;

	for (int i = 0; i < n && correcto; ++i)
	{
		for (int j = 0; j < n && correcto; ++j)
		{
			correcto = c[i] == a[i] + b[i];
		}
	} 

	return correcto;

}
