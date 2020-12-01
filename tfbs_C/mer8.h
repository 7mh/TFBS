#ifndef MER8_H_
#define MER8_H_

typedef struct MER8 {
  char *seq;
  char * seq_compl;
  double e;
  double med;
  double z;
} mer8_t;


void print_8mer(mer8_t m) {
  printf("%s\t%s\t%.2lf\t%.2lf\t%.2lf\n",
	 m.seq, m.seq_compl, m.e, m.med, m.z);
}


void line_to_mer8(char *line, mer8_t *m) {
    // a line is like this (with tabs)
    //8-mer	8-mer	E-score	Median	Z-score
    //AAAAAAAA	TTTTTTTT	0.32779	4276.12	2.0578

    // sequence
    char *s = strtok(line,"\t\n");
    m->seq = s;
    
    // complementary sequence
    s = strtok(NULL,"\t\n");
    m->seq_compl = s;

    // E-score
    s = strtok(NULL,"\t\n");
    m->e = atof(s);

    // median
    s = strtok(NULL,"\t\n");
    m->med = atof(s);

    // Z-score
    s = strtok(NULL,"\t\n");
    m->z = atof(s);
}

#endif