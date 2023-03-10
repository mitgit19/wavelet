
#ifndef F3_H_
#define F3_H_

#include <stdio.h>
#include <stdlib.h>

#include "prng.h"

#include "definitions.h"

#define f3_mul(a, b) (((a) * (b))%3)
#define f3_add(a, b) (((a) + (b))%3)

#define f3_neg(a) (((-(a))+3)%3)

/*
 * a <- a + b
 */

void f3_print(const f3 *a);
void f3_compact_print(const wave_word r0, const wave_word r1);

uint8_t f3_rand(prng_t *PRNG);

uint8_t f3_randnonzero(prng_t *PRNG);

wave_word f3_array_weight(uint8_t *v, int length);

#endif /* F3_H_ */
