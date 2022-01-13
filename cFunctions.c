#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include "myProto.h"
#include <time.h>
//the Conservative Groups
const char * conservative[9] = { "NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
//the Semi-Conservative Groups
const char * semi_conservative[11] = { "SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};

//perform the required task
mutant check(int offset, char * seq1, char * seq2, double * W, int min_max, int start, int end){
   char temp_seq2 [SEQ2_SIZE];
   strcpy(temp_seq2, seq2);
   mutant m;
   double max_score = -100000000, min_score = 100000000, current_score = 0;
   char * p1 = seq1 + offset;//adding the current offset to seq1 
   int index = start;//the required start index for seq2 
   int group;
   char original_c;

   while(seq2[index] != '\0' && index <= end){//continue running the loop until rech end of seq2 or until reaching the required end index for seq2
       char c = 'A';
       original_c = seq2[index];//save the the original char from seq2 for reconstruction
       for (int i = 0 ; i < 26 ; i++){//run through all letters for mutant check
           group = check_group(original_c, c);//check if the current letter (c) can be a mutant according to the Substitution Rules (space or point)
           if(group == POINT || group == SPACE || group == STAR){
               temp_seq2[index] = c;//replace the original letter in a mutant
               current_score = alignment_score(p1, temp_seq2, W);//calculate the alignment score value
               if (min_max){//if searching for maximum score value
                   if(current_score > max_score){
                       max_score = current_score;
                       strcpy(m.seq, temp_seq2);
                       m.score = max_score;
                       m.offset = offset;
		       m.pos = index;
                       m.c = c;
                   }
               }
               else if (current_score < min_score){//if searching for minimum score value
                       min_score = current_score;
                       strcpy(m.seq, temp_seq2);
                       m.score = min_score;
                       m.offset = offset;
		       m.pos = index;
                       m.c = c;
               }
           }
	   strcpy(temp_seq2, seq2);
           c++;
       }
       strcpy(temp_seq2, seq2);
       index++;
   }
   return m;
}

double alignment_score(char * seq1, char * seq2, double * W){//calculating Sequences Alignment Score
    int numStars = 0, numcolons = 0, numPoints = 0, numSpace = 0;
    char * p1 = seq1;
    char * p2 = seq2;
    while(*p2 != '\0' && *p1 != '\0'){
        switch (check_group(*p1, *p2))//check the sign of each pair
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
        p1++;
        p2++;
    }
    return W[0]*numStars - W[1]*numcolons - W[2]*numPoints - W[3]*numSpace;//formula for the Alignment Score
}

int check_group(char c1, char c2){//check the sign of a given pair of letters
    if (c1 == c2)
        return STAR;
    if(check_conservative(c1, c2))
        return COLON;
    else if (check_semi_conservative(c1, c2))
        return POINT;
    else
        return SPACE;
}

int check_conservative(char c1, char c2){//check if a pair of letters are in the same conservative group
    for (int i = 0; i < 9 ; i++)
        if(strchr(conservative[i], c1) != NULL && strchr(conservative[i], c2) != NULL)
            return 1;
    return 0;
}

int check_semi_conservative(char c1, char c2){//check if a pair of letters are in the same semi conservative group
    for (int i = 0; i < 11 ; i++)
        if(strchr(semi_conservative[i], c1) != NULL && strchr(semi_conservative[i], c2) != NULL)
            return 1;
    return 0;
}

void find_optimal(mutant * mutants, double *best_score, mutant *best_mutant, int min_max, int max_offset){//find the mutant with the optimal alignment score from all the best alignment score mutants found
    if (min_max){//if maximum score need to be found
    *best_score = -100000000;
        for (int i = 0;   i < 2*max_offset; i++){
            if(mutants[i].score > *best_score){
                *best_score = mutants[i].score;
                *best_mutant = mutants[i];
            }
		
        }
    }
    else{//if minimum score need to be found
        *best_score = 100000000;
        for (int i = 0;   i < 2*max_offset; i++){
            if(mutants[i].score < *best_score){
                *best_score = mutants[i].score;
                *best_mutant = mutants[i];
            }
        }
    }
} 
