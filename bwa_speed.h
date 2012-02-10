#ifndef BWA_SPEED_H
#define BWA_SPEED_H


/*! 
  @details  Index performance testing.
  */

/*! 
  structure to store the command line options for 'tmap indexspeed'
  */
typedef struct {
    char *fn_fasta;  /*!< the fasta file name (-f) */
    int32_t enum_max_hits;  /*!< the maximum number of hits to enumerate (-e) */
    int32_t kmer_length;  /*!< the kmer length to simulate (-K) */ 
    int32_t kmer_num;  /*!< the number of kmers to simulate (-N) */
    double rand_frac;  /*!< the fraction of random kmers (-R) */
} bwa_speed_opt_t;

/*! 
  main-like function for 'tmap indexspeed'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int bwa_speed(int argc, char *argv[]);

#endif
