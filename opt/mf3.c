#include "mf3.h"
#include <time.h>

mf3* mf3_new(size_t nr_rows, size_t nr_cols) {

	size_t i;
	mf3 *M = calloc(1, sizeof(mf3));
	M->n_rows = nr_rows;
	M->n_cols = nr_cols;
	M->row = (f3_vector*) malloc(M->n_rows * sizeof(f3_vector));
	for (i = 0; i < M->n_rows; ++i) {
		M->row[i] = f3_vector_new(M->n_cols);
	}
	return M;
}

mf3* mf3_transpose(mf3 *M) {
	size_t i;
	mf3 *T = calloc(1, sizeof(mf3));
	T->n_rows = M->n_cols;
	T->n_cols = M->n_rows;
	T->row = (f3_vector*) malloc(T->n_rows * sizeof(f3_vector));
	for (i = 0; i < T->n_rows; i++) {
		T->row[i] = f3_vector_new(T->n_cols);
	}

	//mf3_print(T);

	for (i = 0; i < M->n_rows; i++) {
		for (size_t j = 0; j < M->n_cols; j++) {
			uint8_t a = mf3_coeff(M, i, j);
			mf3_setcoeff(T, j, i, a);
			//mf3_print(T);
		}
	}
	return T;

}

void mf3_free(mf3 *M) {
	for (size_t i = 0; i < M->n_rows; ++i) {
		f3_vector_free(&M->row[i]);
	}
	free(M->row);
	free(M);
}

mf3* mf3_rand(size_t nr_rows, size_t nr_cols, prng_t *PRNG) {

	size_t i;
	mf3 *M = calloc(1, sizeof(mf3));
	M->n_rows = nr_rows;
	M->n_cols = nr_cols;
	M->row = (f3_vector*) malloc(M->n_rows * sizeof(f3_vector));
	for (i = 0; i < M->n_rows; ++i) {
		M->row[i] = f3_vector_new(M->n_cols);
		for (size_t j = 0; j < nr_cols; ++j) {
			f3_vector_setcoeff(&M->row[i], j, rnd_trit(PRNG));
		}
	}

	return M;
}

int mf3_gauss_elim_single(mf3 *M, size_t r, size_t j) {
	size_t i;
	uint8_t a;

	for (i = r; i < M->n_rows; i++) {
		a = f3_vector_coeff(&M->row[i], j);
		if (a != 0) {
			if (a == 2) { // normalize row
				// exchanging r0 and r1 multiplies the row by -1
				f3_vector_neg(&M->row[i]);
				/*wave_word *temp = M->row[i].r0;
				 M->row[i].r0 = M->row[i].r1;
				 M->row[i].r1 = temp;*/

			}

			f3_vector temp = M->row[i];
			M->row[i] = M->row[r];
			M->row[r] = temp;

			for (i = 0; i < M->n_rows; i++) {
				if (i != r) {
					a = f3_vector_coeff(&M->row[i], j);
					if (a == 1)
						f3_vector_sub_inplace(&M->row[i], &M->row[r]);
					else if (a == 2)
						f3_vector_sum_inplace(&M->row[i], &M->row[r]);
				}
			}

			return 1;
		}
	}
	return 0;
}

/***
 support[] is a permutation of {0,..,n-1}.

 The function returns the total number of pivots tried, successful
 or not.

 After the call support[] will have the r=rank(M) ordered
 successful pivots positions first. The last n-r entries of
 support[] contain the n-r other positions.
 ***/
int mf3_gauss_elim(mf3 *M, unsigned int *support) {
	size_t j, k, n, r, d;
	n = M->n_cols;
	k = M->n_rows;
	int pivots[N] = { 0 }; //(int*) malloc(k * sizeof(int));
	int nonpivots[N] = { 0 }; //(int*) malloc(n * sizeof(int));
	for (r = 0, d = 0, j = 0; (j < n) && (r < k); ++j) {
		if (mf3_gauss_elim_single(M, r, support[j]) == 0) {
			nonpivots[d] = support[j];
			d++;
		} else {
			pivots[r] = support[j];
			r++;
		}
	}

	memcpy(support, pivots, r * sizeof(int));
	memcpy(support + r, nonpivots, d * sizeof(int));

	return j;
}

// assumes x.alloc >= y.alloc
void element_copy(f3_vector *x, f3_vector *y) {
	memcpy(x->r0, y->r0, y->alloc * sizeof(wave_word));
	memcpy(x->r1, y->r1, y->alloc * sizeof(wave_word));
}

mf3* mf3_augment(mf3 *H, uint8_t *s) {
	mf3 *Hp = mf3_new(H->n_rows, H->n_cols + 1);
	size_t i = H->n_cols / WORD_LENGTH;
	wave_word z = ((wave_word) 1) << (H->n_cols % WORD_LENGTH);
	for (size_t l = 0; l < H->n_rows; ++l) {
		element_copy(&Hp->row[l], &H->row[l]);
		if (s[l] == 1){
			Hp->row[l].r0[i] ^= z;
		}
		else if (s[l] == 2){
			Hp->row[l].r1[i] ^= z;
			Hp->row[l].r0[i] ^= z;
		}
	}
	// pretend that the last column isn't there
	Hp->n_cols = H->n_cols; // instead of H.coldim + 1
	return Hp;
}

f3_vector f3_vector_copy(f3_vector *y) {
	f3_vector x = f3_vector_new(y->size);
	memcpy(x.r0, y->r0, y->alloc * sizeof(wave_word));
	memcpy(x.r1, y->r1, y->alloc * sizeof(wave_word));

	return x;
}

mf3* mf3_copy(mf3 *M) {

	mf3 *H = calloc(1, sizeof(mf3));
	H->n_rows = M->n_rows;
	H->n_cols = M->n_cols;
	H->row = (f3_vector*) malloc(M->n_rows * sizeof(f3_vector));
	for (size_t i = 0; i < M->n_rows; i++) {
		H->row[i] = f3_vector_copy(&M->row[i]);
	}
	return H;
}

void mf3_mv_mul_v(const mf3 *M, const f3_vector *v, f3_vector *a) {
	size_t i;

	for (i = 0; i < M->n_rows; i++) {
		uint8_t r = f3_vector_mul_sum(&M->row[i], v);
		f3_vector_setcoeff(a, i, r);

	}
}

void mf3_mv_mul(const mf3 *M, const f3_vector *v, uint8_t *a) {
	size_t i;
	for (i = 0; i < M->n_rows; i++) {
		a[i] = f3_vector_dotproduct(&M->row[i], v);
	}
}

void mf3_ma_mul(const mf3 *M, const uint8_t *a, uint8_t *res) {
	f3_vector v = f3_vector_new_copy_from_array(M->row[0].size, a);
	mf3_mv_mul(M, &v, res);
	f3_vector_free(&v);
}

void mf3_print(mf3 *M) {
	size_t i;

	for (i = 0; i < M->n_rows; ++i) {
		printf("\t[");
		f3_vector_print_big(&M->row[i]);
		printf("],\t");
		printf("\n");
	}
	printf("\n");
}


