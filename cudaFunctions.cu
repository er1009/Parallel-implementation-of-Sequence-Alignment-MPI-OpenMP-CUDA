#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "myProto.h"

__global__  void preform_task(mutant * mutants, char * seq1, char * seq2, double * W, int min_max, int start, int end, int maxOffset, char * temp_seq2);
__device__ mutant check_cuda(int offset, char * seq1, char * seq2, double * W, int min_max, int start, int end, int seq2_length, char * temp_seq2);
__device__ double alignment_score_cuda(char * seq1, char * seq2, double * W, char c, int index);
__device__ int check_group_cuda(char c1, char c2);
__device__ int check_conservative_cuda(char c1, char c2);
__device__ int check_semi_conservative_cuda(char c1, char c2);
__device__ int my_strchr(const char *str, char ch);
__device__ char * my_strcpy(char *dest, const char *src);

__global__  void preform_task(mutant * mutants, char * seq1, char * seq2, double * W, int min_max, int start, int end, int maxOffset, int seq2_length, char * temp_seq2) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    // Increment the proper value of the arrray according to thread ID 
    if (i < maxOffset){
        mutants[i] = check_cuda(i, seq1, seq2, W, min_max, start, end, seq2_length, temp_seq2);
    }
}

int computeOnGPU(mutant * mutants, char * seq1, char * seq2, double * W, int min_max, int start, int end, int maxOffset, int seq1_length, int seq2_length){
    char temp_str [SEQ2_SIZE];
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;
    cudaError_t err1 = cudaSuccess;
    cudaError_t err2 = cudaSuccess;
    cudaError_t err3 = cudaSuccess;
    cudaError_t err4 = cudaSuccess;
    cudaError_t err5 = cudaSuccess;

    size_t size_seq1 = (seq1_length + 1) * sizeof(char);
    size_t size_seq2 = (seq2_length + 1) * sizeof(char);
    size_t size_W = 4*sizeof(double);
    size_t size_mutants = maxOffset * sizeof(mutant);

    // Allocate memory on GPU to copy the data from the host
    char *d_seq1;
    char *d_seq2;
    char *d_temp_seq2;
    double * d_W;
    mutant *d_mutants;
    err1 = cudaMalloc((void **)&d_seq1, size_seq1);
    err2 = cudaMalloc((void **)&d_seq2, size_seq2);
    err3 = cudaMalloc((void **)&d_mutants, size_mutants);
    err4 = cudaMalloc((void **)&d_W, size_W);
    err5 = cudaMalloc((void **)&d_temp_seq2, SEQ2_SIZE * sizeof(char));
   

    if (err1 != cudaSuccess || err2 != cudaSuccess || err3 != cudaSuccess || err4 != cudaSuccess ||err5 != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy data from host to the GPU memory
    err1 = cudaMemcpy(d_seq1, seq1, size_seq1, cudaMemcpyHostToDevice);
    err2 = cudaMemcpy(d_seq2, seq2, size_seq2, cudaMemcpyHostToDevice);
    err3 = cudaMemcpy(d_mutants, mutants, size_mutants, cudaMemcpyHostToDevice);
    err4 = cudaMemcpy(d_W, W, size_W, cudaMemcpyHostToDevice);
    err5 = cudaMemcpy(d_temp_seq2, temp_str, SEQ2_SIZE * sizeof(char), cudaMemcpyHostToDevice);

    if (err1 != cudaSuccess || err2 != cudaSuccess || err3 != cudaSuccess || err4 != cudaSuccess || err5 != cudaSuccess) {
        fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Launch the Kernel
    int threadsPerBlock = 256;
    int blocksPerGrid =(maxOffset + threadsPerBlock - 1) / threadsPerBlock;
	
    preform_task<<<blocksPerGrid, threadsPerBlock>>>(d_mutants, d_seq1, d_seq2, d_W, min_max, start, end, maxOffset, seq2_length, d_temp_seq2);
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to launch vectorAdd kernel -  %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy the  result from GPU to the host memory.
    err = cudaMemcpy(mutants, d_mutants, size_mutants, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Free allocated memory on GPU
    if (cudaFree(d_mutants) != cudaSuccess || cudaFree(d_seq1) != cudaSuccess || cudaFree(d_W) != cudaSuccess || cudaFree(d_seq2) != cudaSuccess || cudaFree(d_temp_seq2) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
	
    //create mutant sequence	
	strcpy(temp_str, seq2);
    for(int i = 0 ; i < maxOffset ; i++){
        temp_str[mutants[i].pos] = mutants[i].c;
        strcpy(mutants[i].seq, temp_str);
	strcpy(temp_str, seq2);
    }

    return 0;
}

//perform the required task
__device__ mutant check_cuda(int offset, char * seq1, char * seq2, double * W, int min_max, int start, int end, int seq2_length, char * temp_seq2){
   my_strcpy(temp_seq2, seq2);
   mutant m;
   double max_score = -100000000;
   double min_score = 100000000;
   double current_score = 0;
   char * p1 = seq1 + offset;//adding the current offset to seq1 
   int index = start;//the required start index for seq2
   int group;
   char original_c;
   char c;
   while(seq2[index] != '\0' && index <= end){//continue running the loop until rech end of seq2 or until reaching the required end index for seq2
       c = 'A';
       original_c = seq2[index];//save the the original char from seq2 for reconstruction
       for (int i = 0 ; i < 26 ; i++){//run through all letters for mutant check
           group = check_group_cuda(original_c, c);//check if the current letter (c) can be a mutant according to the Substitution Rules (space or point)
           if(group == POINT ||group == SPACE || group == STAR){
               current_score = alignment_score_cuda(p1, temp_seq2, W, c, index);//calculate the alignment score
               if (min_max){//if searching for maximum alignment score value
                   if(current_score > max_score){
                       max_score = current_score;
                       m.score = max_score;
                       m.offset = offset;
		       m.pos = index;
                       m.c = c;
                   }
               }
               else if (current_score < min_score){//if searching for minimum alignment score value
                       min_score = current_score;
                       m.score = min_score;
                       m.offset = offset;
		       m.pos = index;
                       m.c = c;
               }
           }
           c++;
           
       }
       index++;
   }

   return m;
}

__device__ double alignment_score_cuda(char * seq1, char * seq2, double * W, char c, int index){//calculating Sequences alignment score
    int numStars = 0, numcolons = 0, numPoints = 0, numSpace = 0;
    char * p1 = seq1;
    char * p2 = seq2;
    int pos = 0 ;
    while(*p2 != '\0' && *p1 != '\0'){
	if (pos == index){
		switch (check_group_cuda(*p1, c))//check the sign of each pair
		{
		case STAR:
		    numStars++;
		    break;
		case COLON:
		    numcolons++;
		    break;
		case POINT:
		    numPoints++;
		    break;
		case SPACE:
		    numSpace++;
		    break;
		default:
		    break;
		}
	}
	else{
	        switch (check_group_cuda(*p1, *p2))//check the sign of each pair
		{
		case STAR:
		    numStars++;
		    break;
		case COLON:
		    numcolons++;
		    break;
		case POINT:
		    numPoints++;
		    break;
		case SPACE:
		    numSpace++;
		    break;
		default:
		    break;
		}
	}
        p1++;
        p2++;
	pos++;
    }
    return W[0]*numStars - W[1]*numcolons - W[2]*numPoints - W[3]*numSpace;//formula for the Alignment Score
}

__device__ int check_group_cuda(char c1, char c2){//check the sign of a given pair of latters
    if (c1 == c2)
        return STAR;
    if(check_conservative_cuda(c1, c2))
        return COLON;
    else if (check_semi_conservative_cuda(c1, c2))
        return POINT;
    else
        return SPACE;
}

__device__ int check_conservative_cuda(char c1, char c2){//check if a pair of latters are in the same conservative group
    const char * conservative[9] = { "NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
    for (int i = 0; i < 9 ; i++)
        if(my_strchr(conservative[i], c1) != NULL && my_strchr(conservative[i], c2)!= NULL)
            return 1;
    return 0;
}

__device__ int check_semi_conservative_cuda(char c1, char c2){//check if a pair of latters are in the same semi conservative group
    const char * semi_conservative[11] = { "SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};
    for (int i = 0; i < 11 ; i++)
        if(my_strchr(semi_conservative[i], c1) != NULL  && my_strchr(semi_conservative[i], c2) != NULL)
            return 1;
    return 0;
}

__device__ int my_strchr(const char *str, char ch)//function to preform strchr on cuda
{
for (;; str++) {
        if (*str == ch) return 1;
        if (!*str) return 0;
        }
}

__device__ char * my_strcpy(char *dest, const char *src){
  int i = 0;
  do {
    dest[i] = src[i];}
  while (src[i++] != 0);
  return dest;
}
