/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef BWA_RAND_H
#define BWA_RAND_H

#include <stdint.h>

#define BWA_RAND_NN 312
#define BWA_RAND_MM 156
#define BWA_RAND_UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define BWA_RAND_LM 0x7FFFFFFFULL /* Least significant 31 bits */

#define bwa_rand_drand48(_kr) ((bwa_rand_int(_kr) >> 11) * (1.0/9007199254740992.0))
#define bwa_rand_sample(_kr, _k, _cnt) ((*(_cnt))++ < (_k)? *(_cnt) - 1 : bwa_rand_rand(_kr) % *(_cnt))

/*!
  a thread-safe random number generator
  */
typedef struct {
    int64_t mti;
    uint64_t mt[BWA_RAND_NN];
} bwa_rand_t;

/*!
  @param  seed  the random seed
  @return       the initialized random number generator
 */
bwa_rand_t*
bwa_rand_init(uint64_t seed);

/*!
  @param  r  the initialized random number generator
  @return    a random 64-bit integer
 */
uint64_t 
bwa_rand_int(bwa_rand_t *r);

/*!
  @param  r  the initialized random number generator
  @return    a random double between 0 and 1.
 */
double
bwa_rand_get(bwa_rand_t *r);

/*!
  @param  r  the initialized random number generator
*/
void
bwa_rand_destroy(bwa_rand_t *r);

#endif
