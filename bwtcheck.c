#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"

void bwa_check_core(const char *prefix, int32_t length, int32_t print_sa)
{
  bwt_t *bwt;
  ubyte_t *seq = NULL;
  int32_t i, n;
  uint64_t hash_j;
  int64_t sum;
  int64_t j;
  bwtint_t sa_begin, sa_end;
  char *str;

  str = (char*)calloc(strlen(prefix) + 10, 1);
  strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
  free(str);
  
  /*
  seq = calloc(2, sizeof(ubyte_t));
  seq[0] = 0;
  seq[1] = 1;
  n = bwt_match_exact(bwt, 2, seq, &sa_begin, &sa_end); 
  int k;
  for(k=0;k<2;k++) {
      fputc("ACGTN"[seq[k]], stderr);
  }
  if(0 < n && sa_begin <= sa_end) {
      fprintf(stderr, " %llu %llu %d\n", sa_begin, sa_end, n);
  }
  else {
      fprintf(stderr, " NA NA NA\n");
  }
  */

  for(i=1;i<=length;i++) {
      seq = calloc(i, sizeof(ubyte_t));

      j = 0;
      hash_j = sum = 0;
      while(1) {
          if(i == j) {
              n = bwt_match_exact(bwt, i, seq, &sa_begin, &sa_end); 
              if(0 < n && sa_begin <= sa_end) {
                  sum += n;
              }
              if(1 == print_sa) {
                  int k;
                  for(k=0;k<i;k++) {
                      fputc("ACGTN"[seq[k]], stderr);
                  }
                  if(0 < n && sa_begin <= sa_end) {
                      fprintf(stderr, " %llu %llu %d\n", sa_begin, sa_end, n);
                  }
                  else {
                      fprintf(stderr, " NA NA NA\n");
                  }
              }

              j--;
              while(0 <= j && 3 == seq[j]) {
                  seq[j] = 0;
                  hash_j >>= 2;
                  j--;
              }
              if(j < 0) break;
              seq[j]++;
              hash_j++;
              j++;
          }
          else {
              hash_j <<= 2;
              j++;
          }
      }

      free(seq);

      fprintf(stderr, "kmer:%d sum=%llu expected sum:%llu difference:%lld matched:%s\n",
              i,
              sum,
              bwt->seq_len - i + 1,
              sum - (bwt->seq_len - i + 1),
              sum == (bwt->seq_len - i + 1) ? "true" : "false");
  }

  bwt_destroy(bwt);
}

int bwa_check(int argc, char *argv[])
{
	int c, length = 12, print_sa = 0;

	while ((c = getopt(argc, argv, "l:p")) >= 0) {
		switch (c) {
                  case 'l': length = atoi(optarg); break;
                  case 'p': print_sa = 1; break;
                  default: return 1;
		}
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa check [options] <prefix>\n\n");
		fprintf(stderr, "Options: -l INT    the kmer length to check\n");
		fprintf(stderr, "Options: -p        print out the SA intervals for each kmer\n");
		fprintf(stderr, "\n");
		return 1;
	}
	bwa_check_core(argv[optind], length, print_sa);
	return 0;
}
