#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include "myProto.h"
#include <time.h>

int main(int argc, char *argv[]) {

    int process_num, rank, max_offset, seq1_length = 0, seq2_length = 0, start, end, min_max, process1_best_offset;
    double W[4];
    char seq1 [SEQ1_SIZE];
    char seq2 [SEQ2_SIZE];
    char min_max_str[MINMAX_SIZE];
    double best_score, process1_best_score;
    mutant best_mutant;
    char process1_best_mutant [SEQ2_SIZE];
    mutant * mutants;
    MPI_Status  status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &process_num);
    if (process_num != 2) {
        printf("Run the example with two processes only\n");
        MPI_Abort(MPI_COMM_WORLD, __LINE__);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
   // Divide the tasks between both processes
   if (rank == 0) {
      //read data from file
      FILE* file = fopen(argv[1], "r");
      fscanf(file, "%lf", &W[0]);
      fscanf(file, "%lf", &W[1]);
      fscanf(file, "%lf", &W[2]);
      fscanf(file, "%lf\n", &W[3]);
      fscanf(file, "%s\n", seq1);
      fscanf(file, "%s\n", seq2);
      fscanf(file, "%s\n", min_max_str);
      fclose(file);

      while(seq1[seq1_length] != '\0')//count seq1 length
          seq1_length++;
      
      while(seq2[seq2_length] != '\0')//count seq2 length
          seq2_length++;
      
      max_offset = seq1_length - seq2_length + 1;//calculate maximum offset value
      
      if(strcmp(min_max_str, "maximum") == 0){//check what needded to be found - minimum or maximum value
          min_max = 1;
      }
      else{
          min_max = 0;
      }

      //send data to the second process
      MPI_Send(W, 4, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
      MPI_Send(&seq1_length, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
      MPI_Send(&seq2_length, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
      MPI_Send(seq1, SEQ1_SIZE, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
      MPI_Send(seq2, SEQ2_SIZE, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
      MPI_Send(&max_offset, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
      MPI_Send(&min_max, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);

      //the start and end values of the part in seq2 that need to be handled by the first process
      start = 0;
      end = seq2_length/2;
      mutants = (mutant*)calloc(2*max_offset,sizeof(mutant));//an array of all the best scores mutants that was found for each offset
   } else {
      //recieve data from the first process
      MPI_Recv(W, 4, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&seq1_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&seq2_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(seq1, SEQ1_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(seq2, SEQ2_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&max_offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&min_max, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

      //the start and end values of the part in seq2 that need to be handled by the second process
      start = seq2_length/2;
      end = seq2_length;

      mutants = (mutant*)calloc(2*max_offset,sizeof(mutant));//an array of structs that conain all best alignment scores mutants that was found for each offset
}

   // On each process - perform a first half of its task with OpenMP     
#pragma omp parallel for
   for (int i = 0;   i < max_offset; i++)     
      mutants[i] = check(i, seq1, seq2, W, min_max, start, (start + end)/2);       


   // On each process - perform a second half of its task with CUDA
   if (computeOnGPU(mutants + max_offset, seq1, seq2, W, min_max, (start + end)/2, end, max_offset, seq1_length, seq2_length) != 0)
    MPI_Abort(MPI_COMM_WORLD, __LINE__);
   
   // On each process - find the mutant with the optimal alignment score from all the best alignment score mutants found
   find_optimal(mutants, &best_score, &best_mutant, min_max, max_offset);

   // Collect the result on one of processes
   if (rank == 0){
      //Collect the results from the second process
      MPI_Recv(&process1_best_score, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(process1_best_mutant, SEQ2_SIZE, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&process1_best_offset, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
      
      FILE* result_file = fopen("./result.txt", "w");
      //check for the best score and write to file
      if((min_max && best_mutant.score > process1_best_score) || (!min_max && best_mutant.score < process1_best_score)){
            fprintf(result_file, "mutant: %s\n", best_mutant.seq);
            fprintf(result_file, "offset: %d score: %lf", best_mutant.offset, best_mutant.score);
      }
      else{
            fprintf(result_file, "mutant: %s\n", process1_best_mutant);
            fprintf(result_file, "offset: %d score: %lf", process1_best_offset, process1_best_score);
	
      }

      fclose(result_file);
      printf("finish\n");fflush(stdout);
   } 
   else{
      //send the results to the first process
      MPI_Send(&best_mutant.score, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(best_mutant.seq, SEQ2_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&best_mutant.offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); 
   } 

   MPI_Finalize();
   return 0;
}
