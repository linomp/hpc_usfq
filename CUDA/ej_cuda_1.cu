#include <stdio.h>
#include <stdlib.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//The number of character in the encrypted text
#define N 1024

void checkCUDAError(const char*);
int get_text_length(const char * fname);
void read_file(const char*, int*);
void write_file(const char*, int*, int );

#define A 15
#define B 27
#define M 128
#define A_MMI_M 111

__device__ int modulo(int a, int b){
	int r = a % b;
	r = (r < 0) ? r + b : r;
	return r;
}

__global__ void decrypt_kernel(int *d_input, int *d_output, int length)
{
	for (int i = 0; i< length; i++){
		char x = d_input[i];
		d_output[i] = 111*(x-B) % M; 
	}
}

__global__ void decrypt_multiblock_kernel(int *d_input, int *d_output, int length)
{
	// ...
}


int main(int argc, char *argv[])
{
	int *h_input, *h_output;
	int *d_input, *d_output;
	unsigned int size;

	const char * fname;

	if (argc < 2) printf("Debe ingresar el nombre del archivo\n");
	else
		fname = argv[argc-1];

	int length = get_text_length(fname);

	size = length * sizeof(int);

	// reservo memoria para h_input y h_output
	h_input = (int *)malloc(size);
	h_output = (int *)malloc(size);

	// reservar memoria en la GPU para d_input y d_output
	cudaMalloc(&d_input, size); cudaMalloc(&d_output, size);
	checkCUDAError("Memory allocation");

	// leo el archivo con el mensaje cifrado
	read_file(fname, h_input);

	// transferir el arreglo de entrada al dispositivo
	cudaMemcpy(d_input, h_input, size, cudaMemcpyHostToDevice);
	checkCUDAError("Input transfer to device");

	// configurar la grilla de threads
	dim3 blocksPerGrid (1, 1, 1) 
	dim3 threadsPerBlock (N, 1, 1)

	// ejecutar el kernel
	decrypt_kernel <<< blocksPerGrid, threadsPerBlock >>>( d_input, d_output, length );

	cudaThreadSynchronize();
	checkCUDAError("Kernel execution");


	// transferir el contenido de d_output a la memoria de la CPU
	cudaMemcpy(h_output, d_output, size, cudaMemcpyDeviceToHost);
	checkCUDAError("Result transfer to host");

	printf("Texto desencriptado:\n");

	for (int i = 0; i < length; i++) {
		printf("%c", (char)h_output[i]); 
	}
	printf("\n");

	// liberar memoria en el dispositivo para d_input y d_output
	cudaFree(d_input); cudaFree(d_output);
	checkCUDAError("Free memory");

	/* free host buffers */
	free(h_input);
	free(h_output);

	return 0;
}


void checkCUDAError(const char *msg)
{
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "CUDA ERROR: %s: %s.\n", msg, cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
}


	
int get_text_length(const char * fname)
{
	FILE *f = NULL;
	f = fopen(fname, "r"); 

	size_t pos = ftell(f);    
	fseek(f, 0, SEEK_END);    
	size_t length = ftell(f); 
	fseek(f, pos, SEEK_SET);  

	fclose(f);

	return length;
}

void read_file(const char * fname, int* input)
{
	// printf("leyendo archivo %s\n", fname );

	FILE *f = NULL;
	f = fopen(fname, "r"); 
	if (f == NULL){
		fprintf(stderr, "Error: Could not find %s file \n", fname);
		exit(1);
	}

	int c; 
	while ((c = getc(f)) != EOF) {
		*(input++) = c;
	}

	fclose(f);
}

void write_file(const char * fname, int* input, int length)
{
	FILE *f = NULL;
	f = fopen(fname, "w"); 
	if (f == NULL){
		fprintf(stderr, "Error: Could not find %s file \n", fname);
		exit(1);
	}

	for (int i = 0; i < length; ++i)
	{
		putc((char)input[i],f);
	}

	fclose(f);
}