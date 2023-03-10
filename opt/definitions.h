#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <stdio.h>
#include <stdint.h>
#include "fips202.h"

#include "params.h"

typedef uint64_t wave_word;

#define WORD_LENGTH (8*sizeof(wave_word))
#define HASH_SIZE 64

typedef struct {
	wave_word r0;
	wave_word r1;
} f3;

typedef struct {
	wave_word *r0;
	wave_word *r1;
	size_t alloc, size;
} f3_vector;

typedef struct {
	f3_vector *row;
	size_t n_rows, n_cols;
} mf3;

typedef struct prng {
	uint8_t buff[256];
	uint8_t buff2;
	uint8_t buff3;
	uint8_t available;
	uint8_t available2;
	uint8_t available3;
	unsigned long bytecount;
	uint8_t init_str[25];
	shake256incctx sha3;
} prng_t;

typedef struct {
	uint8_t a[N2];
	uint8_t b[N2];
	uint8_t c[N2];
	uint8_t d[N2];
	uint8_t x[N2];
} coeff_t;

typedef struct {
	mf3 *HU;
	mf3 *HV;
	mf3 *Sinv;
	coeff_t coeff;
	unsigned int perm[N];
} wave_sk_t;

typedef mf3 *wave_pk_t;

#endif /* DEFINITIONS_H_ */
