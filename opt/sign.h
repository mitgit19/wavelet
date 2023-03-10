
#ifndef SRC_WAVE_SIGN_H_
#define SRC_WAVE_SIGN_H_

#include <stdlib.h>
#include <string.h>

#include "fips202.h"

#include "f3.h"
#include "vf3.h"
#include "mf3.h"

#include "prng.h"
#include "randombytes.h"

#include "definitions.h"

#include "commons.h"

#include "wave.h"
#include "params.h"

#include "hash.h"

// v a vector of size 2 * length
// number of non-zero pairs (v[i], v[i+length])
int nz_pairs(uint8_t *v, int length);

// v is a vector of size length.
// support will contain all integers from 0 to length-1.
// Non-zero postions of v come first in support.
// Side effect: weight is set to the Hamming weight of v.
void supp(uint8_t *v, int length, int *weight, unsigned int *suppV);

void decodeV(uint8_t *eV, uint8_t *sV, mf3 *HV, prng_t *PRNG);
//int8_t* decodeV(int8_t *sV, mf3_t HV, prng_t PRNG);

void decodeU(uint8_t *e, uint8_t *sU, uint8_t *eV, mf3 *HU, wave_sk_t *sk,
		prng_t *PRNG);
//int8_t* decodeU(int8_t *sU, int8_t *eV, mf3_t HU, wave_sk_t sk, prng_t PRNG);

//void sign_wave(int8_t* e, int8_t *s, wave_sk_t sk, prng_t PRNG);

void sign(f3_vector *signature, f3_vector *message_hash, uint8_t *salt,
		const uint8_t *message, const size_t mlen, wave_sk_t *sk);

#endif /* SRC_WAVE_SIGN_H_ */
