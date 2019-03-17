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
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
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

	int length = sizeof(v1)/sizeof(v1[0]);
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
  cudaMemcpy(d_input_1, v1, size, cudaMemcpyHostToDevice);   
	cudaMemcpy(d_input_2, v2, size, cudaMemcpyHostToDevice); 

  // configurar la grilla de threads
	dim3 blocksPerGrid ( (int) ceil(length/N), 1, 1) ;
	dim3 threadsPerBlock (N, 1, 1);

	// ejecutar el kernel
	add_elements_kernel <<< blocksPerGrid, threadsPerBlock >>>( d_input_1, d_input_2, d_output, length );

	// sólo para medir tiempos, porque el memcoy ya sincroniza internamente
	cudaThreadSynchronize(); 

	// transferir el contenido de d_output a la memoria de la CPU
	cudaMemcpy(h_output, d_output, size, cudaMemcpyDeviceToHost); 

	printf("Texto desencriptado:\n");

	for (int i = 0; i < length; i++) {
		printf("%.2f", (float)h_output[i]); 
	}
	printf("\n");

	// liberar memoria en el dispositivo para d_input y d_output
  cudaFree(d_input_1); 
  cudaFree(d_input_2);
  cudaFree(d_output); 
  
	/* free host buffers */
	//free(v1);
	//free(v2);

	return 0;
} 