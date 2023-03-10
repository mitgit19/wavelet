#ifndef MF3_H_
#define MF3_H_

#include "definitions.h"

#include "vf3.h"

#include "prng.h"

mf3* mf3_new(size_t nr_rows, size_t nr_cols);

void mf3_free(mf3 *M);

mf3* mf3_rand(size_t nr_rows, size_t nr_cols, prng_t *PRNG);

static inline void mf3_setcoeff(mf3 *M, int i, int j, int8_t a) {
	f3_vector_setcoeff(&M->row[i], j, a);
}

static inline uint8_t mf3_coeff(mf3 *M, int i, int j) {
	return f3_vector_coeff(&M->row[i], j);
}

void mf3_times_vector(const mf3 *M, const f3_vector *v, f3_vector *res);

int mf3_gauss_elim(mf3 *M, unsigned int *support);

void mf3_print(mf3 *M);

mf3* mf3_copy(mf3 *M);

mf3* mf3_augment(mf3 *H, uint8_t *s);

void mf3_ma_mul(const mf3 *M, const uint8_t *a, uint8_t *res);

void mf3_mv_mul_v(const mf3 *M, const f3_vector *v, f3_vector *a);

void mf3_mv_mul(const mf3 *M, const f3_vector *v, uint8_t *a);

mf3* mf3_transpose(mf3 *M);

#endif /* MF3_H_ */
