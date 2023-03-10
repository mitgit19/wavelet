#ifndef VF3_H_
#define VF3_H_

#define popcount __builtin_popcountll

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <immintrin.h>

#include "libpopcnt.h"

#include "definitions.h"

#include "prng.h"

f3_vector f3_vector_new(size_t length);

f3_vector f3_vector_new_copy_from_array(const size_t length, const uint8_t *a);

size_t f3_vector_weight(const f3_vector *vec);

void f3_vector_free(f3_vector *x);

void f3_vector_zero(f3_vector *x);

void f3_vector_setcoeff(f3_vector *x, size_t j, uint8_t a);

uint8_t f3_vector_coeff_iszero(f3_vector *x, int j);

uint8_t f3_vector_coeff_isone(f3_vector *x, int j);

uint8_t f3_vector_coeff_istwo(f3_vector *x, int j);

uint8_t f3_vector_get_coeff(f3_vector *x, int j);

uint8_t f3_vector_mul_sum(const f3_vector *x, const f3_vector *y);

void f3_vector_print(const f3_vector *a);

void f3_vector_print_big(const f3_vector *a);

uint8_t f3_vector_sum(f3_vector *x);

uint8_t f3_vector_dotproduct(const f3_vector *x, const f3_vector *y);

void f3_vector_mul(const f3_vector *x, const f3_vector *y, f3_vector *res);

void f3_vector_rand(f3_vector *a, prng_t *PRNG);

void f3_vector_set_from_array(f3_vector *x, const uint8_t *v,
		const size_t length);

static inline void f3_vector_sum_inplace(f3_vector *a, const f3_vector *b) {
	wave_word tmp = 0;
	for (size_t i = 0; i < a->alloc; i++) {
		tmp = 0;
		tmp = (a->r0[i] ^ b->r1[i]) & (a->r1[i] ^ b->r0[i]);
		a->r0[i] = (a->r0[i] ^ b->r0[i]) | ((a->r1[i] ^ b->r0[i]) ^ b->r1[i]);
		a->r1[i] = tmp;
	}
}

static inline void f3_vector_sum_inplace_f(wave_word *a_r0, wave_word *a_r1,
		wave_word *b_r0, wave_word *b_r1, size_t alloc) {
	wave_word tmp = 0;
	for (size_t i = 0; i < alloc; i++) {
		tmp = 0;
		tmp = (a_r0[i] ^ b_r1[i]) & (a_r1[i] ^ b_r0[i]);
		a_r0[i] = (a_r0[i] ^ b_r0[i]) | ((a_r1[i] ^ b_r0[i]) ^ b_r1[i]);
		a_r1[i] = tmp;
	}
}

static inline void f3_vector_add(const f3_vector *a, const f3_vector *b,
		f3_vector *c) {
	for (size_t i = 0; i < a->alloc; i++) {
		wave_word tmp = (a->r1[i] ^ b->r0[i]);
		c->r1[i] = (a->r0[i] ^ b->r1[i]) & tmp;
		c->r0[i] = (a->r0[i] ^ b->r0[i]) | (tmp ^ b->r1[i]);
	}
}
static inline void f3_vector_neg_vector(const f3_vector *a, f3_vector *b) {
	memcpy(b->r0, a->r0, a->alloc * sizeof(wave_word));
	memcpy(b->r1, a->r1, a->alloc * sizeof(wave_word));
	for (int i = 0; i < b->alloc; i++) {
		b->r1[i] = b->r1[i] ^ b->r0[i];
	}

}
static inline void f3_vector_neg(f3_vector *b) {
	for (int i = 0; i < b->alloc; i++) {
		b->r1[i] = b->r1[i] ^ b->r0[i];
	}
}

static inline void f3_vector_copy_vector(const f3_vector *a, f3_vector *b) {
	memcpy(b->r0, a->r0, a->alloc * sizeof(wave_word));
	memcpy(b->r1, a->r1, a->alloc * sizeof(wave_word));

}
// c = a - b
static inline void f3_vector_sub(const f3_vector *a, const f3_vector *b,
		f3_vector *c) {

	for (size_t i = 0; i < a->alloc; i++) {
		wave_word tmp = a->r0[i] ^ b->r0[i];
		c->r1[i] = (tmp ^ b->r1[i]) & (a->r1[i] ^ b->r0[i]);
		c->r0[i] = tmp | ((a->r1[i] ^ b->r1[i]));
	}
}

static inline void f3_vector_sub_inplace(f3_vector *a, const f3_vector *b) {

	for (size_t i = 0; i < a->alloc; i++) {
		wave_word tmp_t = (a->r0[i] ^ (b->r0[i] ^ b->r1[i]))
				& (a->r1[i] ^ b->r0[i]);
		a->r0[i] = (a->r0[i] ^ b->r0[i]) | ((a->r1[i] ^ b->r1[i]));
		a->r1[i] = tmp_t;
	}

}

static inline void f3_vector_sub_inplace_f(wave_word *a_r0, wave_word *a_r1,
		wave_word *b_r0, wave_word *b_r1, size_t alloc) {

	for (size_t i = 0; i < alloc; i++) {
		wave_word tmp_t = (a_r0[i] ^ (b_r0[i] ^ b_r1[i])) & (a_r1[i] ^ b_r0[i]);
		a_r0[i] = (a_r0[i] ^ b_r0[i]) | ((a_r1[i] ^ b_r1[i]));
		a_r1[i] = tmp_t;
	}

}

static inline uint8_t f3_vector_coeff(const f3_vector *x, int j) {
	int i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	if (((x->r0[i] >> j) & 1) && ((x->r1[i] >> j) & 1)) {
		return 2;
	} else if (((x->r0[i] >> j) & 1)) {
		return 1;
	} else {
		return 0;
	}

}

static inline void f3_vector_sum_inplace_avx(const f3_vector *a,
		const f3_vector *b) {

	int offset = 0;
	for (offset = 0; offset < a->alloc - 4; offset = offset + 4) {
		double *addr_a = (double*) (&a->r0[offset]);
		__m256d a_l = _mm256_loadu_pd(addr_a);
		__m256d a_h = _mm256_loadu_pd((double*) a->r1 + offset);
		__m256d b_l = _mm256_loadu_pd((double*) b->r0 + offset);
		__m256d b_h = _mm256_loadu_pd((double*) b->r1 + offset);

		__m256d res_1 = _mm256_xor_pd(a_l, b_h);
		__m256d res_2 = _mm256_xor_pd(a_h, b_l);
		__m256d c_high = _mm256_and_pd(res_1, res_2);
		_mm256_storeu_pd((double*) a->r1 + offset, c_high);

		res_1 = _mm256_xor_pd(a_l, b_l);
		res_2 = _mm256_xor_pd(a_h, b_l);
		__m256d res_3 = _mm256_xor_pd(res_2, b_h);
		__m256d c_low = _mm256_or_pd(res_1, res_3);
		//__m256i* ptr_l = (__m256i*) (c->r0 + ptr);
		_mm256_storeu_pd((double*) a->r0 + offset, c_low);

	}
	if (offset < a->alloc) {
		if (a->alloc % offset == 2) {
			double *addr_a = (double*) (&a->r0[offset]);
			__m128d a_l = _mm_load_pd(addr_a);
			__m128d a_h = _mm_load_pd((double*) a->r1 + offset);
			__m128d b_l = _mm_load_pd((double*) b->r0 + offset);
			__m128d b_h = _mm_load_pd((double*) b->r1 + offset);

			__m128d res_1 = _mm_xor_pd(a_l, b_h);
			__m128d res_2 = _mm_xor_pd(a_h, b_l);
			__m128d c_high = _mm_and_pd(res_1, res_2);
			_mm_storeu_pd((double*) a->r1 + offset, c_high);

			res_1 = _mm_xor_pd(a_l, b_l);
			res_2 = _mm_xor_pd(a_h, b_l);
			__m128d res_3 = _mm_xor_pd(res_2, b_h);
			__m128d c_low = _mm_or_pd(res_1, res_3);
			//__m256i* ptr_l = (__m256i*) (c->r0 + ptr);
			_mm_storeu_pd((double*) a->r0 + offset, c_low);
		} else if (a->alloc % offset == 1) {
			wave_word tmp = 0;
			tmp = (a->r0[offset] ^ b->r1[offset])
					& (a->r1[offset] ^ b->r0[offset]);
			a->r0[offset] = (a->r0[offset] ^ b->r0[offset])
					| ((a->r1[offset] ^ b->r0[offset]) ^ b->r1[offset]);
			a->r1[offset] = tmp;
		}
	}
}

static inline void f3_vector_sub_inplace_avx(f3_vector *a, const f3_vector *b) {
	int offset = 0;
	for (offset = 0; offset < a->alloc - 4; offset = offset + 4) {
		__m256d a_l = _mm256_loadu_pd((double*) a->r0 + offset);
		__m256d a_h = _mm256_loadu_pd((double*) a->r1 + offset);

		__m256d b_l = _mm256_loadu_pd((double*) b->r0 + offset);
		__m256d b_h = _mm256_loadu_pd((double*) b->r1 + offset);

		__m256d res_1 = _mm256_xor_pd(a_l, b_l);

		__m256d res_2 = _mm256_xor_pd(a_h, b_h);

		__m256d c_low = _mm256_or_pd(res_1, res_2);

		_mm256_storeu_pd((double*) a->r0 + offset, c_low);

		res_1 = _mm256_xor_pd(b_l, b_h);

		__m256d res_4 = _mm256_xor_pd(a_h, b_l);

		res_2 = _mm256_xor_pd(a_l, res_1);

		__m256d c_high = _mm256_and_pd(res_4, res_2);

		_mm256_storeu_pd((double*) a->r1 + offset, c_high);

	}

	if (offset < a->alloc) {
		if (a->alloc - offset == 2) {
			__m128d a_l = _mm_load_pd((double*) a->r0 + offset);
			__m128d a_h = _mm_load_pd((double*) a->r1 + offset);
			__m128d b_l = _mm_load_pd((double*) b->r0 + offset);
			__m128d b_h = _mm_load_pd((double*) b->r1 + offset);

			__m128d res_1 = _mm_xor_pd(a_l, b_l);

			__m128d res_2 = _mm_xor_pd(a_h, b_h);

			__m128d c_low = _mm_or_pd(res_1, res_2);

			_mm_storeu_pd((double*) a->r0 + offset, c_low);

			res_1 = _mm_xor_pd(b_l, b_h);

			__m128d res_4 = _mm_xor_pd(a_h, b_l);

			res_2 = _mm_xor_pd(a_l, res_1);

			__m128d c_high = _mm_and_pd(res_4, res_2);

			_mm_storeu_pd((double*) a->r1 + offset, c_high);

		}
	}
}

static inline void f3_vector_sum_inplace_avx_f(wave_word *a_r0, wave_word *a_r1,
		wave_word *b_r0, wave_word *b_r1, size_t alloc) {

	int offset = 0;
	for (offset = 0; offset < alloc - 4; offset = offset + 4) {
		double *addr_a = (double*) (&a_r0[offset]);
		__m256d a_l = _mm256_loadu_pd(addr_a);
		__m256d a_h = _mm256_loadu_pd((double*) a_r1 + offset);
		__m256d b_l = _mm256_loadu_pd((double*) b_r0 + offset);
		__m256d b_h = _mm256_loadu_pd((double*) b_r1 + offset);

		__m256d res_1 = _mm256_xor_pd(a_l, b_h);
		__m256d res_2 = _mm256_xor_pd(a_h, b_l);
		__m256d c_high = _mm256_and_pd(res_1, res_2);
		_mm256_storeu_pd((double*) a_r1 + offset, c_high);

		res_1 = _mm256_xor_pd(a_l, b_l);
		res_2 = _mm256_xor_pd(a_h, b_l);
		__m256d res_3 = _mm256_xor_pd(res_2, b_h);
		__m256d c_low = _mm256_or_pd(res_1, res_3);
		//__m256i* ptr_l = (__m256i*) (c->r0 + ptr);
		_mm256_storeu_pd((double*) a_r0 + offset, c_low);

	}
	if (offset < alloc) {
		if (alloc % offset == 2) {
			double *addr_a = (double*) (&a_r0[offset]);
			__m128d a_l = _mm_load_pd(addr_a);
			__m128d a_h = _mm_load_pd((double*) a_r1 + offset);
			__m128d b_l = _mm_load_pd((double*) b_r0 + offset);
			__m128d b_h = _mm_load_pd((double*) b_r1 + offset);

			__m128d res_1 = _mm_xor_pd(a_l, b_h);
			__m128d res_2 = _mm_xor_pd(a_h, b_l);
			__m128d c_high = _mm_and_pd(res_1, res_2);
			_mm_storeu_pd((double*) a_r1 + offset, c_high);

			res_1 = _mm_xor_pd(a_l, b_l);
			res_2 = _mm_xor_pd(a_h, b_l);
			__m128d res_3 = _mm_xor_pd(res_2, b_h);
			__m128d c_low = _mm_or_pd(res_1, res_3);
			//__m256i* ptr_l = (__m256i*) (c->r0 + ptr);
			_mm_storeu_pd((double*) a_r0 + offset, c_low);
		} else if (alloc % offset == 1) {
			wave_word tmp = 0;
			tmp = (a_r0[offset] ^ b_r1[offset]) & (a_r1[offset] ^ b_r0[offset]);
			a_r0[offset] = (a_r0[offset] ^ b_r0[offset])
					| ((a_r1[offset] ^ b_r0[offset]) ^ b_r1[offset]);
			a_r1[offset] = tmp;
		}
	}
}

static inline void f3_vector_sub_inplace_avx_f(wave_word *a_r0, wave_word *a_r1,
		wave_word *b_r0, wave_word *b_r1, size_t alloc) {
	int offset = 0;
	for (offset = 0; offset < alloc - 4; offset = offset + 4) {
		__m256d a_l = _mm256_loadu_pd((double*) a_r0 + offset);
		__m256d a_h = _mm256_loadu_pd((double*) a_r1 + offset);

		__m256d b_l = _mm256_loadu_pd((double*) b_r0 + offset);
		__m256d b_h = _mm256_loadu_pd((double*) b_r1 + offset);

		__m256d res_1 = _mm256_xor_pd(a_l, b_l);

		__m256d res_2 = _mm256_xor_pd(a_h, b_h);

		__m256d c_low = _mm256_or_pd(res_1, res_2);

		_mm256_storeu_pd((double*) a_r0 + offset, c_low);

		res_1 = _mm256_xor_pd(b_l, b_h);

		__m256d res_4 = _mm256_xor_pd(a_h, b_l);

		res_2 = _mm256_xor_pd(a_l, res_1);

		__m256d c_high = _mm256_and_pd(res_4, res_2);

		_mm256_storeu_pd((double*) a_r1 + offset, c_high);

	}

	if (offset < alloc) {
		if (alloc - offset == 2) {
			__m128d a_l = _mm_load_pd((double*) a_r0 + offset);
			__m128d a_h = _mm_load_pd((double*) a_r1 + offset);
			__m128d b_l = _mm_load_pd((double*) b_r0 + offset);
			__m128d b_h = _mm_load_pd((double*) b_r1 + offset);

			__m128d res_1 = _mm_xor_pd(a_l, b_l);

			__m128d res_2 = _mm_xor_pd(a_h, b_h);

			__m128d c_low = _mm_or_pd(res_1, res_2);

			_mm_storeu_pd((double*) a_r0 + offset, c_low);

			res_1 = _mm_xor_pd(b_l, b_h);

			__m128d res_4 = _mm_xor_pd(a_h, b_l);

			res_2 = _mm_xor_pd(a_l, res_1);

			__m128d c_high = _mm_and_pd(res_4, res_2);

			_mm_storeu_pd((double*) a_r1 + offset, c_high);

		}
	}
}

f3_vector f3_vector_one_new(size_t length);

#endif /* VF3_H_ */
