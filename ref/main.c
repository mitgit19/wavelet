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
//size_t m4r_opt_k(size_t n) {
  //size_t k = log2(n);
  //return k;
//}

// Power 2^k (Considering that for now the implementation is for F2)
size_t two_pow(size_t k) {
   size_t res = 2 << (k-1);
  return res;
}

// MakeTable implementation
void make_table(mf3 *A, size_t r_str, size_t c_str, size_t k, mf3 *T){

	// TODO: Replace the hardcoded section with separate function
	// Using the pre-processed "fundamental indexes" (the indexes correspond to T | start from 0) 
	size_t I[k];
	memset(I, 0, k * sizeof(size_t));

	I[0] = 1;
	I[1] = 2;
	I[2] = 4;
	I[3] = 8;
	I[4] = 16;
	I[5] = 32;
	I[6] = 64;
	I[7] = 128;
	I[8] = 256;
	I[9] = 512;
	I[10] = 1024;
	I[11] = 2048;

	// Setting the first two rows of T after the 0th one
	T->row[1] = A->row[r_str + (k - 1)];

	// Added to ensure that we can work with 1 pivot
	if (k > 1) {
		T->row[2] = A->row[r_str + (k - 2)];
	}

	// TODO: Replace the hardcoded section with separate function
	// Indicator used for keeping track with the upcoming fundamental index
	size_t i = 2;

	// Storing the last fundamental row (to avoid multiple indexing)
	f3_vector fund_row = T->row[2];

	// Storing the number of rows in the prior to the fundamental that will be useful ([0..(I[i]-1)] rule)
	size_t dep_i = 1; // For I[i] = 2 that will be 1

	// Creating T based on 2^k additions (considering min k = 1)
	for (size_t j = 3; j < two_pow(k); j++) {

		// If we are on a fundamental index -> increment by one the indicator
		if (j == I[i]) {
			T->row[j] = A->row[r_str + (k - i - 1)];
			fund_row = T->row[j];
			dep_i = j-1;

			// Ensuring that we won't go out of bounds
			if ((k - 1) > i) {
				i++;
			}
		}

		else {
			// j increments by one but dep_i remains constant thus we effectively start from the row that is dep_row rows above the last fundamental row and move downwards until the row prior to the fundamental itself
			f2_vector_add(&T->row[j - dep_i - 1],&fund_row,&T->row[j]);
		}
	}
}

// AddRowsFromTable implementation
void add_rows_from_table(mf3 *A, size_t r_str, size_t r_end, size_t c_str, size_t k, mf3 *T) {
	
	// Within the set region of rows
	for (size_t r = r_str; r < r_end; r++) {
		size_t index = 0;

		// For the first k columns starting from c_str
		for (size_t j = 0; j < k; j++) {
			index = index + (f3_vector_get_coeff(&A->row[r], c_str + j) << (k - j - 1));
		}

		// Adding the row of T corresponding to the previous index value to the row r of A
		f2_vector_sum_inplace(&A->row[r], &T->row[index]);
	}
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
	//size_t k = m4r_opt_k(m);
	
	size_t k = 6;

	while(c < n) {

		// If out of boundary (nu of columns)
		if (c+k > n) {
			k = n-c;
		}

		// Computing GaussSubmatrix
		size_t k_bar = f2_m4r_mf3_gauss_elim(A, k, r, c);

		// TODO: Remove - Used for debugging
		// printf("Pivots: %d\n", k_bar);
		// printf("Position r: %d\n", r);
		// printf("Position c: %d\n", c);

		// TODO: Remove - Used for debugging
		// printf("Eliminated A: \n");
		// mf3_print(A);

		// If pivots found
		if (k_bar > 0) {

			// Computing 2^k
			size_t two_pow_k = two_pow(k_bar);

			// Creating table T (size 2^k x n)
			mf3* T = mf3_new(two_pow_k, n);

			// Creating integer array L (size 2^k)
		//	size_t L[two_pow_k];
		//	memset(L, 0, two_pow_k * sizeof(size_t));

			// Making the table with all the linear combs for 2^k
			make_table(A,r,c,k_bar,T);

			// TODO: Remove - Used for debugging
			// printf("Creating the table T: \n");
			// mf3_print(T);
			
			// Adding the rows from the table
			add_rows_from_table(A,0,r,c,k_bar,T);
			add_rows_from_table(A,r+k_bar,m,c,k_bar,T);
		}

		// TODO: Remove - Used for debugging
		// printf("Added rows A: \n");
		// mf3_print(A);


		// Updating r,c
		r += k_bar;
		c += k_bar;

		// if (k != k_bar) { // I believe this step is wrongly typed in the M4R algorithm (I saw the GitHub repo has also removed it)
		// 	c += 1;
		// }

		if (k_bar == 0){
			c +=1;
		}
	}
}


#include "debug.h"

int main(void) {

	init(PARAMS_ID);

	printf("------\t Testing region for M4R.\t--------\n");

	// Generating random matrices
	prng_t *PRNG;
	PRNG = prng_init(2);
	//const size_t mat_nu_rows = 2887;
	//const size_t mat_nu_cols = 8492;
	 const size_t mat_nu_rows = 13;
	 const size_t mat_nu_cols = 17;

	mf3* H1 = f2_mf3_rand(mat_nu_rows, mat_nu_cols, PRNG);
	mf3* H2 = mf3_copy(H1);

	// Visualizing the random matrix
	// mf3_print(H1);

	// printf("%d\n", H1->n_rows);
	// printf("%d\n", H1->n_cols);

	// Timing the original Gaussian Elimination

	// Creating the support vector
	unsigned int support[17] = { 0 };
	for (int i = 0; i < mat_nu_cols; ++i)
		support[i] = i;

	clock_t t1;
    t1 = clock();

	size_t piv_tried = f2_mf3_gauss_elim(H1,support);

	t1 = clock() - t1;
	double time_taken = ((double)t1)/CLOCKS_PER_SEC;
	printf("The original Gaussian Elim took %f seconds to execute \n", time_taken);

	// Visualizing the result obtained by the original gaussian elimination
	//mf3_print(H1);

	// Timing the M4R Gaussian Elimination
	clock_t t2;
    t2 = clock();

	m4r_impl(H2);
	t2 = clock() - t2;
	time_taken = ((double)t2)/CLOCKS_PER_SEC;
	printf("The M4R Gaussian Elim took %f seconds to execute \n", time_taken);

	// Visualizing the result obtained by the M4R gaussian elimination
	//mf3_print(H2);

	// Testing the makeTable function
	// mf3* T = mf3_new(two_pow(3), 10);
	// make_table(H1,0,0,3,T);
	// mf3_print(T);

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
