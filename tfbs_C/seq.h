#ifndef SEQ_H_
#define SEQ_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char *read_seq(const char *filename) {
  FILE * f = fopen(filename, "r");
  assert(f != NULL);
  
  long length;
  fseek (f, 0, SEEK_END);
  length = ftell (f);
  fseek (f, 0, SEEK_SET);

  char *buffer = (char *) malloc (length + 1);
  assert(buffer != NULL);
  fread (buffer, 1, length, f);
  fclose (f);

  // now go through and condense, get rid of whitespace
  int curr=0, to=0;
  for(; curr < length; ) {
    if (isspace(buffer[curr])) {
      curr++;
    }
    else {
      buffer[to] = buffer[curr];
      curr++;
      to++;
    }
  }
  buffer[to] = '\0';
  buffer = (char *) realloc(buffer, to+1);

  return buffer;
}

// convert sequence to integer, viewing the
// sequence as a base-4 "number" with A=0, C=1, G=2, T=3
typedef enum {SEQ_A=0, SEQ_C, SEQ_G, SEQ_T} DNA_BP;
unsigned seq_to_i(char *s, int bps) {
  if (s == NULL) return -1;
  
  unsigned total = 0;
  for(; *s && bps > 0; s++) {
    total <<= 2;
    unsigned bb = 0;
    switch(*s) {
    case 'a':
    case 'A':
      bb = SEQ_A;
      break;
    case 't':
    case 'T':
    case 'u':
    case 'U':
      bb = SEQ_T;
      break;
    case 'c':
    case 'C':
      bb = SEQ_C;
      break;
    case 'g':
    case 'G':
      bb = SEQ_G;
      break;
    }
    total |= bb;
    bps--; // how many basepairs to keep reading
  }
  if (bps > 0) // if didn't get to 8 bp's before end of string
    return -1; 
  return total;
}

// return the code of the reverse complement of seq
// note that the code for acgt were chosen above
// so that the complement of the code is the code for
// the complementary bp
unsigned reverse_compl(unsigned seq, int bps) {
  unsigned result = 0;
  for(int i=0; i < bps; i++) {
    result <<= 2; // 2 bits per bp
    result |= ((seq & 0x3) ^ 0x3); // compl to get compl base
    seq >>= 2;
  }
  return result;
}


// return the code of mutating seq by changing the bp at
// pos to be val.  pos is 0 on the right of the sequence
unsigned mutate(unsigned seq, int pos, unsigned val) {
  seq &= (~ (0x3 << (2*pos)));
  seq |= (val << (2*pos));
  return seq;
}

#endif