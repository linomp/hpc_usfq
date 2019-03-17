#include <stdio.h>
#include <stdlib.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//The number of character in the encrypted text
#define N 1024

void checkCUDAError(const char*);
int get_text_length(const char * fname); 
 
__global__ void add_elements_kernel(float *d_input_1, float *d_input_2, float *d_output, int length)
{
	// qué thread soy?
	// aplicando offset porque hay más bloques
  //int idx = blockDim.x * blockIdx.x + threadIdx.x;
  int idx = threadIdx.x;
	// verificar, xq pueden lanzarse más hilos que elementos
	if( idx < length ){ 
		float res = d_input_1[idx] + d_input_2[idx];
		// calcular y guardar valor en indice correspondiente en el buffer de salida
		d_output[idx] = res; 
	}
}


int main(int argc, char *argv[])
{
  float v1[5] = {1000.0, 2.0, 3.4, 7.0, 50.0};
  float v2[5] = {1.0, 4.0, 4, 1, 50.0};
  float *h_output;
	float *d_input_1, *d_input_2, *d_output;
	unsigned int size; 
  
	float t_deviceToHost, t_kernel, t_hostToDevice;
  struct timeval t_i, t_dth, t_k, t_htd; 

  int length = sizeof(v1)/sizeof(float); 

	size = length * sizeof(float);

	// reservo memoria para h_input y h_output
	//v1 = (float *)malloc(size);
  //v2 = (float *)malloc(size);
  h_output = (float *)malloc(size);

	// reservar memoria en la GPU para d_input y d_output
  cudaMalloc(&d_input_1, size); 
  cudaMalloc(&d_input_2, size); 
  cudaMalloc(&d_output, size); 

  // transferir el arreglo de entrada al dispositivo
  gettimeofday(&t_i, NULL);
  cudaMemcpy(d_input_1, v1, size, cudaMemcpyHostToDevice);   
  cudaMemcpy(d_input_2, v2, size, cudaMemcpyHostToDevice); 
  gettimeofday(&t_htd, NULL);

  t_hostToDevice = (double)t_htd.tv_sec * 1000.0 + (double)t_htd.tv_usec / 1000.0 -
           ((double)t_i.tv_sec * 1000.0 + (double)t_i.tv_usec / 1000.0);

  // configurar la grilla de threads
  //dim3 blocksPerGrid ( (int) ceil(length/N), 1, 1) ;
  dim3 blocksPerGrid (1, 1, 1) ;
	dim3 threadsPerBlock (N, 1, 1);

  // ejecutar el kernel
  gettimeofday(&t_i, NULL);
	add_elements_kernel <<< blocksPerGrid, threadsPerBlock >>>( d_input_1, d_input_2, d_output, length );

	// sólo para medir tiempos, porque el memcopy ya sincroniza internamente
  cudaThreadSynchronize(); 
  gettimeofday(&t_k, NULL);
  t_kernel = (double)t_k.tv_sec * 1000.0 + (double)t_k.tv_usec / 1000.0 -
           ((double)t_i.tv_sec * 1000.0 + (double)t_i.tv_usec / 1000.0);

	// transferir el contenido de d_output a la memoria de la CPU
  gettimeofday(&t_i, NULL);
  cudaMemcpy(h_output, d_output, size, cudaMemcpyDeviceToHost); 
  gettimeofday(&t_dth, NULL);
  t_deviceToHost = (double)t_dth.tv_sec * 1000.0 + (double)t_dth.tv_usec / 1000.0 -
           ((double)t_i.tv_sec * 1000.0 + (double)t_i.tv_usec / 1000.0);

	printf("Suma (GPU):\n");
	for (int i = 0; i < length; i++) {
		printf("%.2f, ", (float)h_output[i]); 
	}
  printf("\n");
  
  printf("Suma (Host):\n");
  verifyInHost(v1, v2, length);

  printf("\nTiempo transf. host-to-device %f ms\nTiempo kernel %f ms\nTiempo transf. device-to-host % f ms\n ", t_hostToDevice, t_kernel, t_deviceToHost);


	// liberar memoria en el dispositivo para d_input y d_output
  cudaFree(d_input_1); 
  cudaFree(d_input_2);
  cudaFree(d_output); 
  
	/* free host buffers */
	//free(v1);
	//free(v2);

	return 0;
} 

void verifyInHost(float* v1, float* v2, int length){
  for (int i = 0; i < length; i++) {
		printf("%.2f, ", v1[i] + v2[i]); 
  }
  printf("\n");
}