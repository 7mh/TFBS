#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#include "mer8.h"
#include "seq.h"

#define NUM_8MERS 65536
double e[NUM_8MERS];


void mutate(int ind[], int c, char * seq);

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("  Usage: ./tf8merAlign SEQ_FILE\n");
    exit(0);
  }

  // read the sequence
  char *seq = read_seq(argv[1]);
  //printf("%s", seq);

  // initialize e-score array
  for(int i=0; i < NUM_8MERS; i++)
    e[i] = 0.0;

  // read 8mers lines and put into e array
  mer8_t m; 
  char *line = NULL; size_t n = 0;
  int header = 1;
  while (getline(&line, &n, stdin) > 0) {
    if (header) {
      header = 0;
      continue;
    }
    line_to_mer8(line, &m);

    // convert sequence to index into e array, put e in there
    unsigned seq_code, rev_compl;
    seq_code = seq_to_i(m.seq, 8);
    rev_compl = reverse_compl(seq_code, 8);
    assert(seq_code < NUM_8MERS);
    assert(rev_compl < NUM_8MERS);
    e[seq_code] = m.e;
    e[rev_compl] = m.e;
    
    //printf("  %x, %x | ", seq_code, rev_compl);
    //print_8mer(m);
  }

  // scan through sequence note places of alignment with escore high
  char *seq_targets = strdup(seq);
  int run_length = 0;
  
  int len;
  len = strlen(seq);
  int  ind[len];         //create array of len equal to seq file
  memset(ind, 0, sizeof(ind));

  int c =0 ;

  for(int i=0; seq[i] != '\0'; i++) {
    int seq_code = seq_to_i(&seq[i], 8);
    if (seq_code < 0) continue;
    
    assert(seq_code < NUM_8MERS);
    
    if (e[seq_code] > 0.3) {
      run_length++;
      seq_targets[i] = '0' + run_length; // flag it as a place with a match
    }
    else {
        if (run_length > 0) {
	    // put total run length on all the bp's of run that just ended
	        int j = i-1;
	        while (j >= 0 && seq_targets[j] != ' ') {
	            seq_targets[j] = '0' + run_length;
                //getting bp's indx where match is greater than 4
                if (run_length > 3){                
                    ind[c] = j; c++;
                }
	            j--;
                // insert -1 when match sequence ends
                if (seq_targets[j] == ' '){
                    ind[c] = -1; c++;
                }
	        }

        }
        // reset run length
        run_length = 0;
        seq_targets[i] = ' ';
    }
  }

  printf("%s\n", seq_targets);

  printf("Total matched 8mers = %d\n",len,c);
  for (int i =0; i < c; i++){
      printf("%d,",ind[i]);
  }
  
  printf("\n");
  
  mutate( ind, c, seq);
    


  free(seq);
  free(seq_targets);
  return 0;
}

void mutate(int ind[], int c, char * seq){
    char* t;
    char wrd [4] = {'T','C','G','A'};
    int code_t;

    for (int i =0; i < c; i++){   //iter over all seq (match len) >= 4
      if (ind[i] == -1)
          continue;

      t = strndup(&seq[ind[i]], 8 );
      code_t = seq_to_i(t, 8);
      printf("%s, %d, e=%lf\n", t, code_t, e[code_t]);
      
      int mismatch_found_flg = 0;
      char hold;

      for (int j = 0; j < 8; j++){
          hold = t[j];  
          mismatch_found_flg = 0;
          
          for (int k = 0; k < 4; k++){
                if(t[j] == wrd[k])
                    continue;
                t[j] = wrd[k];
                code_t = seq_to_i(t, 8);
                if(e[code_t] < 0.3){
                    printf("Mismatch found ind= %d seq=%s new e=%lf \n",j,t,e[code_t]);
                    mismatch_found_flg = 1;
                    break;
                }
            }
            t[j] = hold;
            if(mismatch_found_flg)
                break;
      }
    }
}

