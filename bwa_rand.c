#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "bwa_rand.h"

/**
 * 64-bit Mersenne Twister pseudorandom number generator. Adapted from:
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
 * which was written by Takuji Nishimura and Makoto Matsumoto and released
 * under the 3-clause BSD license.
*/

static void 
bwa_rand_srand0(uint64_t seed, bwa_rand_t *r)
{
  r->mt[0] = seed;
  for (r->mti = 1; r->mti < BWA_RAND_NN; ++r->mti)
    r->mt[r->mti] = 6364136223846793005ULL * (r->mt[r->mti - 1] ^ (r->mt[r->mti - 1] >> 62)) + r->mti;
}

bwa_rand_t *
bwa_rand_init(uint64_t seed)
{
  bwa_rand_t *r;
  r = calloc(1, sizeof(bwa_rand_t));
  bwa_rand_srand0(seed, r);
  return r;
}

uint64_t 
bwa_rand_int(bwa_rand_t *r)
{
  uint64_t x;
  static const uint64_t mag01[2] = { 0, 0xB5026F5AA96619E9ULL };
  if (r->mti >= BWA_RAND_NN) {
      int i;
      if (r->mti == BWA_RAND_NN + 1) bwa_rand_srand0(5489ULL, r);
      for (i = 0; i < BWA_RAND_NN - BWA_RAND_MM; ++i) {
          x = (r->mt[i] & BWA_RAND_UM) | (r->mt[i+1] & BWA_RAND_LM);
          r->mt[i] = r->mt[i + BWA_RAND_MM] ^ (x>>1) ^ mag01[(int)(x&1)];
      }
      for (; i < BWA_RAND_NN - 1; ++i) {
          x = (r->mt[i] & BWA_RAND_UM) | (r->mt[i+1] & BWA_RAND_LM);
          r->mt[i] = r->mt[i + (BWA_RAND_MM - BWA_RAND_NN)] ^ (x>>1) ^ mag01[(int)(x&1)];
      }
      x = (r->mt[BWA_RAND_NN - 1] & BWA_RAND_UM) | (r->mt[0] & BWA_RAND_LM);
      r->mt[BWA_RAND_NN - 1] = r->mt[BWA_RAND_MM - 1] ^ (x>>1) ^ mag01[(int)(x&1)];
      r->mti = 0;
  }
  x = r->mt[r->mti++];
  x ^= (x >> 29) & 0x5555555555555555ULL;
  x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
  x ^= (x << 37) & 0xFFF7EEE000000000ULL;
  x ^= (x >> 43);
  return x;
}

double
bwa_rand_get(bwa_rand_t *r)
{
  return bwa_rand_int(r) / (double)UINT64_MAX;
}

void
bwa_rand_destroy(bwa_rand_t *r)
{
  free(r);
}
