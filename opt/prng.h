
#ifndef PRNG_H_
#define PRNG_H_

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "params.h"
#include "fips202.h"

#include "definitions.h"

#include "yates.h"

uint8_t rnd8(prng_t *PRNG);

uint8_t rnd_bit(prng_t *PRNG);

uint8_t rnd_trit(prng_t *PRNG);

uint8_t random_trit_nonzero(prng_t *PRNG);

uint16_t rnd16(prng_t *PRNG);

uint16_t rnd_short(int n, prng_t *PRNG);

prng_t* prng_init(unsigned long seed);

void prng_clear(prng_t *PRNG);

void randperm(unsigned int *a, int n);

void rand_shuffle(unsigned int *a, int n, int t, prng_t *PRNG);

#endif /* PRNG_H_ */
