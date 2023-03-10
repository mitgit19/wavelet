#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "cpucycles.h"
#include "vf3.h"
#include "mf3.h"
#include "prng.h"
#include "commons.h"
#include "compress.h"

#include "keygen.h"
#include "sign.h"
#include "verify.h"

#define MAXMLEN 16

void print_signature(f3_vector *sig) {

	f3_vector_print(sig);
	printf("\n");
}

#define ITERATIONS 1

double timekey_f[ITERATIONS];
double time_sign_f[ITERATIONS];
double time_verify_f[ITERATIONS];

void swap(double *a, double *b) {
	double t = *a;
	*a = *b;
	*b = t;
}

// A function to implement bubble sort
void bubbleSort(double arr[], int n) {
	int i, j;
	for (i = 0; i < n - 1; i++)

		// Last i elements are already in place
		for (j = 0; j < n - i - 1; j++)
			if (arr[j] > arr[j + 1])
				swap(&arr[j], &arr[j + 1]);
}


// (WIP) - Temporary section (TODO: Create separate .c and corresponding .h files)

// Deducing the optimal number for k
size_t m4r_opt_k(size_t n) {
  size_t k = log2(n);
  return k;
}

// Power 2^k (Considering that for now the implementation is for F2)
size_t two_pow(size_t k) {
  size_t res = 2 << k;
  return res;
}

// MakeTable implementation
void make_table(mf3 *A, size_t r_str, size_t c_str, size_t k, mf3 *T, size_t *L) {
	// TODO
}

// AddRowsFromTable implementation
void add_rows_from_table(mf3 *A, size_t r_str, size_t r_end, size_t c_str, size_t k, mf3 *T, size_t *L) {
	// TODO
}

/* Adding M4RI implementation
Input: A - mxn matrix.
Input: k - an integer k>0. For our case this will be deduced automatically.
Result: A in reduced row echelon form
*/

void m4r_impl(mf3 *A){

	// Initialising both r,c to 0
	size_t r = 0;
	size_t c = 0;

	// Finding number of rows and columns
	size_t m = A->n_rows;
	size_t n = A->n_cols;

	// Finding k based on n_columns (since n>m)
	size_t k = m4r_opt_k(m);

	while(c < n) {

		// If out of boundary (nu of columns)
		if (c+k > n) {
			k = n-c;
		}

		// Computing GaussSubmatrix
		size_t k_bar = mf3_gauss_elim_single(A, r, c);

		// If pivots found
		if (k_bar > 0) {

			// Computing 2^k
			size_t two_pow_k = two_pow(k_bar);

			// Creating table T (size 2^k x n)
			mf3* T = mf3_new(two_pow_k, n);

			// Creating integer array L (size 2^k)
			size_t L[two_pow_k];
			memset(L, 0, two_pow_k * sizeof(size_t));

			// Making the table with all the linear combs for 2^k
			make_table(A,r,c,k_bar,T,L);
			
			// Adding the rows from the table
			add_rows_from_table(A,0,r,c,k_bar,T,L);
			add_rows_from_table(A,r+k_bar,m,c,k_bar,T,L);
		}

		// Updating r,c
		r += k_bar;
		c += k_bar;

		if (k != k_bar) {
			c += 1;
		}
	}
}


#include "debug.h"

int main(void) {

	init(PARAMS_ID);

	printf("------\t Warming up.\t--------\n");

		wave_pk_t pk;
		wave_sk_t sk;
		uint8_t salt[SALT_SIZE] = { 0 };

		f3_vector signature = f3_vector_new(N);

		f3_vector m_hash = f3_vector_new(N - K);
		uint8_t mi[MAXMLEN] = { 0 };

		keygen(&sk, &pk);
	//	save_matrix(pk, "pk1.bin");

	//	read_matrix(pk, "pk1.bin");
		f3_vector_zero(&signature);
		f3_vector_zero(&m_hash);
		memset(mi, 0, MAXMLEN * sizeof(uint8_t));
		randombytes(mi, MAXMLEN);

		sign(&signature, &m_hash, salt, mi, MAXMLEN, &sk);

		size_t alloc = (1 + ((K - 1)) / WORD_LENGTH);
		wave_word a_r0[(1 + ((K - 1)) / WORD_LENGTH)] = { 0 };
		wave_word a_r1[(1 + ((K - 1)) / WORD_LENGTH)] = { 0 };

		memcpy(a_r0, signature.r0 + 45, 88 * sizeof(wave_word));
		memcpy(a_r1, signature.r1 + 45, 88 * sizeof(wave_word));
		shift7(a_r0, alloc);
		shift7(a_r1, alloc);

		int verify = Nverify(salt, mi, MAXMLEN, a_r0, a_r1, pk);
		printf("works? %d\n", verify);

		wave_pk_clear(pk);
		wave_sk_clear(&sk);
		f3_vector_free(&signature);
		f3_vector_free(&m_hash);
	printf("\t FINISHED \n");
	cleanup();
	return EXIT_SUCCESS;
}
