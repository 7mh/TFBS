#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "mer8.h"

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("  Usage: ./tf8merFilter eScoreCutoff\n");
    exit(0);
  }
  double e_cutoff = atof(argv[1]);

  mer8_t m;
  
  // read lines
  char *line = NULL; size_t n = 0;
  int header = 1;
  while (getline(&line, &n, stdin) > 0) {
    if (header) {
      header = 0;
      printf("%s", line);
      continue;
    }
    line_to_mer8(line, &m);

    if (m.e >= e_cutoff)
      print_8mer(m);
  }
  
  return 0;
}
