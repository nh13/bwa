/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/shm.h>
#include <time.h>
#include "bwt.h"
#include "bntseq.h"
#include "bwa_rand.h"
#include "bwa_speed.h"

static bwtint_t 
bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, bwtint_t sapos, int len, int *strand)
{
  bwtint_t pos_f;
  int is_rev;
  pos_f = bns_depos(bns, bwt_sa(bwt, sapos), &is_rev); // pos_f
  *strand = !is_rev;
  /* NB: For gapped alignment, pacpos may not be correct, which will be fixed
   * in refine_gapped_core(). This line also determines the way "x" is
   * calculated in refine_gapped_core() when (ext < 0 && is_end == 0). */
  if (is_rev) pos_f = pos_f + 1 < len? 0 : pos_f - len + 1; // mapped to the forward strand
  return pos_f; // FIXME: it is possible that pos_f < bns->anns[ref_id].offset
}


static int32_t  
bwa_speed_test(bwt_t *bwt, bntseq_t *bns, ubyte_t *pacseq, clock_t *total_clock, bwa_speed_opt_t *opt)
{
  bwa_rand_t *rand;
  int32_t i, j, num_found, found;
  bwtint_t k, l;
  bwtint_t pacpos;
  uint8_t *seq = NULL;
  clock_t start_clock = 0;

  rand = bwa_rand_init(13);
  seq = malloc(opt->kmer_length * sizeof(uint8_t));

  num_found = 0;
  for(i=0;i<opt->kmer_num;i++) {
      if(0 < i && 0 == (i % 50000)) {
          fprintf(stderr, "processed %d kmers\n", i);
      }
      if(bwa_rand_get(rand) < opt->rand_frac) {
          for(j=0;j<opt->kmer_length;j++) {
              seq[j] = (uint8_t)(bwa_rand_get(rand) * 4);
          }
      }
      else {
          // get a position (one-based)
          pacpos = 1 + (bwtint_t)(bwa_rand_get(rand) * (bns->l_pac - opt->kmer_length + 1));
          for(j=0;j<opt->kmer_length;j++) {
              seq[j] = pacseq[pacpos>>2] >> ((~pacpos&3)<<1) & 3;
          }
      }
      found = 0;
      start_clock = clock();
      if(0 < bwt_match_exact(bwt, opt->kmer_length, seq, &k, &l)) {
          if(0 <= opt->enum_max_hits && (l - k + 1) <= opt->enum_max_hits) {
              while(k<=l) {
                  int strand;
                  pacpos = bwa_sa2pos(bns, bwt, k, opt->kmer_length, &strand);
                  k++;
              }
          }
          found = 1;
      }
      (*total_clock) += clock() - start_clock;
      if(0 < found) {
          num_found++;
      }
  }
  fprintf(stderr, "processed %d kmers\n", i);

  bwa_rand_destroy(rand);
  free(seq);

  return num_found;
}

static void 
bwa_speed_core(bwa_speed_opt_t *opt)
{
  bwt_t *bwt=NULL;
  bntseq_t *bns=NULL;
  clock_t total_clock= 0;
  time_t start_time, end_time;
  int32_t num_found;
  char *str = NULL;
  ubyte_t *pacseq = NULL;

  // load bwt
  str = (char*)calloc(strlen(opt->fn_fasta) + 10, 1);
  strcpy(str, opt->fn_fasta); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
  strcpy(str, opt->fn_fasta); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
  free(str);
  str = NULL;

  // load bns
  str = (char*)calloc(strlen(opt->fn_fasta) + 10, 1);
  strcpy(str, opt->fn_fasta);
  //strcat(strcpy(str, opt->fn_fasta), ".nt");
  bns = bns_restore(str);
  pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
  rewind(bns->fp_pac);
  fread(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);
  free(str);
  str = NULL;

  // clock on
  start_time = time(NULL);

  // run the speed test
  num_found = bwa_speed_test(bwt, bns, pacseq, &total_clock, opt);

  // clock off
  end_time = time(NULL);

  // free
  bwt_destroy(bwt);
  bns_destroy(bns);
  free(pacseq);

  // print the results
  fprintf(stderr, "[bwt speed] found %d out of %d (%.2f%%)\n", num_found, opt->kmer_num, (100.0 * num_found) / opt->kmer_num);
  fprintf(stderr, "[bwt speed] wallclock time: %d seconds\n", (int)difftime(end_time,start_time));
  fprintf(stderr, "[bwt speed] index lookup cpu cycle time: %.2f seconds\n", (float)(total_clock) / CLOCKS_PER_SEC);
}

static int 
usage(bwa_speed_opt_t *opt)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: bwa speed [options]");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options (required):\n");
  fprintf(stderr, "         -f FILE     the FASTA file name to test [%s]\n", 
                    (NULL == opt->fn_fasta) ? "not using" : opt->fn_fasta);
  fprintf(stderr, "         -e INT      the maximum number of hits to enumerate (-1 for unlimited, 0 to disable) [%d]\n", opt->enum_max_hits);
  fprintf(stderr, "         -K INT      the kmer length to simulate [%d]\n", opt->kmer_length);
  fprintf(stderr, "         -N INT      the number of kmers to simulate [%d]\n", opt->kmer_num);
  fprintf(stderr, "         -R FLOAT    the fraction of random kmers [%.2lf]\n", opt->rand_frac);
  fprintf(stderr, "Options (optional):\n");
  fprintf(stderr, "         -h          print this message\n");
  fprintf(stderr, "\n");
  return 1;
}

int 
bwa_speed(int argc, char *argv[])
{
  int c;
  bwa_speed_opt_t opt;

  opt.fn_fasta = NULL;
  opt.enum_max_hits = 1024;
  opt.kmer_length = 12;
  opt.kmer_num = 100000;
  opt.rand_frac = 0.0;

  while((c = getopt(argc, argv, "f:e:K:N:R:h")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = strdup(optarg); break;
        case 'e':
          opt.enum_max_hits = atoi(optarg); break;
        case 'K':
          opt.kmer_length = atoi(optarg); break;
        case 'N':
          opt.kmer_num = atoi(optarg); break;
        case 'R':
          opt.rand_frac = atof(optarg); break;
        case 'h':
        default:
          return usage(&opt);
      }
  }

  if(argc != optind || 1 == argc) {
      return usage(&opt);
  }
  if(NULL == opt.fn_fasta) {
      fprintf(stderr, "required option -f");
      exit(1);
  }
  if(opt.kmer_length <= 0) {
      fprintf(stderr, "the kmer length to simulate must be greater than zero (-K)");
      exit(1);
  }
  if(opt.kmer_num <= 0) {
      fprintf(stderr, "the numer of kmers to simulate must be greater than zero (-N)");
      exit(1);
  }
  if(opt.rand_frac < 0 || 1 < opt.rand_frac) {
      fprintf(stderr, "the option -R must be between 0 and 1");
      exit(1);
  }

  bwa_speed_core(&opt);

  free(opt.fn_fasta);

  return 0;
}
