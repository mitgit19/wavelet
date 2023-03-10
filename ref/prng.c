#include "prng.h"

void fill(prng_t *PRNG) {
	//sha3_512(PRNG->buff, PRNG->buff, HASH_SIZE);
	shake256_inc_squeeze(PRNG->buff, 64, &PRNG->sha3);
	PRNG->available = 64;
}

void fill2(prng_t *PRNG) {
	PRNG->buff2 = rnd8(PRNG);
	PRNG->available2 = 8;
}

void fill3(prng_t *PRNG) {
	// If the random byte is in [243,255], the first 5 trits of its
	// decomposition in base 3 are not uniformly distributed. We discard
	// those bytes (5%)
	PRNG->buff3 = 254;
	while (PRNG->buff3 >= 243) {
		PRNG->buff3 = rnd8(PRNG);
	}
	; // 243 = 3^5
	PRNG->available3 = 5;
}

uint8_t rnd8(prng_t *PRNG) {
	PRNG->bytecount++;
	if (PRNG->available == 0) {
		fill(PRNG);
	}
	PRNG->available--;
	return PRNG->buff[PRNG->available];

}

uint8_t rnd_bit(prng_t *PRNG) {
	uint8_t c;
	if (PRNG->available2 == 0) {
		fill2(PRNG);
	}
	PRNG->available2--;
	c = PRNG->buff2 & 1;
	PRNG->buff2 >>= 1;
	return c;
}

uint8_t rnd_trit(prng_t *PRNG) {
	uint8_t c;
	if (PRNG->available3 == 0) {
		fill3(PRNG);
	}
	PRNG->available3--;
	c = ((PRNG->buff3 - 1) + 3) % 3;
	PRNG->buff3 /= 3;
	return c;
}

uint8_t random_trit_nonzero(prng_t *PRNG) {
	uint8_t c = 0;
	do {
		if (PRNG->available3 == 0) {
			fill3(PRNG);
		}
		PRNG->available3--;
		c = ((PRNG->buff3 - 1) + 3) % 3;
		PRNG->buff3 /= 3;
	} while (c == 0);
	return c;
}

uint16_t rnd16(prng_t *PRNG) {
	uint16_t r = rnd8(PRNG);
	r <<= 8;
	r ^= rnd8(PRNG);
	return r;
}

/* Random integer uniform between 0 and n-1 included.
 We assume n < 2^16 (else equivalent to rnd16())
 */
uint16_t rnd_short(int n, prng_t *PRNG) {
	uint16_t r, max;
	/* (max + 1) = l * n, the largest multiple of n <= 2^16 */
	max = 0xffff - (0x10000 % n);
	/* accept any r in [0, max] */
	/* if n >= 2^16 then max = -1 = 2^16 - 1 */
	do {
		r = rnd16(PRNG);
	} while (r > max); // else r is not uniformly distributed
	return r % n;
}

prng_t* prng_init(unsigned long seed) {
	prng_t *PRNG = (prng_t*) calloc(1, sizeof(prng_t));

	// 20 is the decimal size of the largest 'unsigned long'
	//PRNG->init_str = malloc(strlen(type) + 4 + 20 + 1);
	memcpy(PRNG->init_str, PARAMS_ID, sizeof(uint8_t) * PARAMS_ID_SIZE);
	memcpy(PRNG->init_str + PARAMS_ID_SIZE, &seed, sizeof(unsigned long) * 1);
	shake256incctx sha3;

	shake256_inc_init(&sha3);

	shake256_inc_absorb(&sha3, PRNG->init_str, 25);
	shake256_inc_squeeze(PRNG->buff, 64, &sha3);
	PRNG->sha3 = sha3;
	shake256_inc_squeeze(PRNG->buff, 128, &PRNG->sha3);

	//PRNG->sha3 = &sha3;
	PRNG->available = HASH_SIZE;
	PRNG->available2 = 0;
	PRNG->available3 = 0;
	PRNG->bytecount = 0;
	PRNG->buff3 = PRNG->buff[0];
	PRNG->buff2 = PRNG->buff[1];
	return PRNG;
}

void prng_clear(prng_t *PRNG) {
	//free(PRNG->init_str);
	//shake256_inc_ctx_release(PRNG->sha3);
	//free(PRNG->sha3->ctx);
	free(PRNG);
}

int* identity_perm(int n) {
	int *a = (int*) malloc(n * sizeof(int));
	for (int i = 0; i < n; i++)
	a[i] = i;
	return a;
}

void randperm(unsigned int *l, int n) {
	for (int i = 0; i < n; i++) {
		l[i] = i;
	}
	fisher_yates_shuffle(l, n);

}

/* randomly shuffle the array a[] of length n so that the first t
 coordinates is a random uniform subset of size t of the array
 values */
void rand_shuffle(unsigned int *a, int n, int t, prng_t *PRNG) {
	int i, r, c;

	for (i = 0; i < t; i++) {
		r = i + rnd_short(n - i, PRNG); /* i <= r < n uniform */
		c = a[r];
		a[r] = a[i];
		a[i] = c;
	}
}
