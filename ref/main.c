#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
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

#define ITERATIONS 40

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

#include "debug.h"

int main(void) {

	init(PARAMS_ID);

	// printf("------\t Testing region for M4R.\t--------\n");

	// Generating random matrices
	// prng_t *PRNG;
	// PRNG = prng_init(2);
	// const size_t mat_nu_rows = 2887;
	// const size_t mat_nu_cols = 8492;
	// // const size_t mat_nu_rows = 13;
	// // const size_t mat_nu_cols = 17;

	// mf3* H1 = mf3_rand(mat_nu_rows, mat_nu_cols, PRNG);
	// mf3* H2 = mf3_copy(H1);

	// Visualizing the random matrix
	// // mf3_print(H1);

	// printf("%d\n", H1->n_rows);
	// printf("%d\n", H1->n_cols);

	// Timing the original Gaussian Elimination

	// Creating the support vector
	// unsigned int support[2887] = { 0 };
	// for (int i = 0; i < mat_nu_rows; ++i)
	// 	support[i] = i;

	// clock_t t1;
    // t1 = clock();

	// size_t piv_found_1 = mf3_gauss_elim(H1,support);

	// t1 = clock() - t1;
	// double time_taken = ((double)t1)/CLOCKS_PER_SEC;
	// printf("The original Gaussian Elim took %f seconds to execute \n", time_taken);

	// Visualizing the result obtained by the original gaussian elimination
	// // mf3_print(H1);

	// Timing the M4R Gaussian Elimination
	// clock_t t2;
    // t2 = clock();

	// m4r_impl(H2);
	// t2 = clock() - t2;
	// time_taken = ((double)t2)/CLOCKS_PER_SEC;
	// printf("The M4R Gaussian Elim took %f seconds to execute \n", time_taken);

	// Visualizing the result obtained by the M4R gaussian elimination
	// // mf3_print(H2);

	// Checking whether the results match
	// printf("The results obtained match: %s \n", (mf3_are_equal(H1,H2)) ? "true" : "false");

		printf("------\t Warming up.\t--------\n");

	printf("\nWithout M4R in keygen: \n");
	clock_t start, end;
	memset(timekey_f, 0, sizeof(double) * ITERATIONS);
	memset(time_sign_f, 0, sizeof(double) * ITERATIONS);
	memset(time_verify_f, 0, sizeof(double) * ITERATIONS);
	
	for (int i = 0; i < ITERATIONS; i++) {
		wave_pk_t pk;
		wave_sk_t sk;
		uint8_t salt[SALT_SIZE] = { 0 };

		f3_vector signature = f3_vector_new(N);

		f3_vector m_hash = f3_vector_new(N - K);
		uint8_t mi[MAXMLEN] = { 0 };

		start = clock();
		keygen(&sk, &pk);
		end = clock();
		timekey_f[i] = (double) (end - start) * 1000.0 / CLOCKS_PER_SEC;
	//	save_matrix(pk, "pk1.bin");

	//	read_matrix(pk, "pk1.bin");
		f3_vector_zero(&signature);
		f3_vector_zero(&m_hash);
		memset(mi, 0, MAXMLEN * sizeof(uint8_t));
		randombytes(mi, MAXMLEN);

		start = clock();
		sign(&signature, &m_hash, salt, mi, MAXMLEN, &sk);
		end = clock();
		time_sign_f[i] = (double) (end - start) * 1000.0 / CLOCKS_PER_SEC;

		size_t alloc = (1 + ((K - 1)) / WORD_LENGTH);
		wave_word a_r0[(1 + ((K - 1)) / WORD_LENGTH)] = { 0 };
		wave_word a_r1[(1 + ((K - 1)) / WORD_LENGTH)] = { 0 };

		memcpy(a_r0, signature.r0 + 45, 88 * sizeof(wave_word));
		memcpy(a_r1, signature.r1 + 45, 88 * sizeof(wave_word));
		shift7(a_r0, alloc);
		shift7(a_r1, alloc);

		start = clock();
		int verify = Nverify(salt, mi, MAXMLEN, a_r0, a_r1, pk);
		end = clock();
		time_verify_f[i] = (double) (end - start) * 1000.0 / CLOCKS_PER_SEC;
		//printf("works? %d\n", verify);

		wave_pk_clear(pk);
		wave_sk_clear(&sk);
		f3_vector_free(&signature);
		f3_vector_free(&m_hash);
	}

	bubbleSort(timekey_f, ITERATIONS - 1);
	bubbleSort(time_verify_f, ITERATIONS - 1);
	bubbleSort(time_sign_f, ITERATIONS - 1);
	double sign_avg = 0;
	double verify_avg = 0;
	double keygen_avg = 0;
	for (int i = 0; i < ITERATIONS; i++) {
		sign_avg += time_sign_f[i];
		verify_avg += time_verify_f[i];
		keygen_avg += timekey_f[i];
	}

	printf("keygen() MEDIAN %f milli-seconds \n", timekey_f[ITERATIONS / 2]);

	printf("keygen() AVG %f milli-seconds \n", (keygen_avg / ITERATIONS));

	printf("sign() MEDIA %f milli-seconds\n", time_sign_f[ITERATIONS / 2]);

	printf("sign() AVG took %f milli-seconds \n", (sign_avg / ITERATIONS));

	printf("Nverify() MEDIAN %f milli-seconds \n",
			time_verify_f[ITERATIONS / 2]);

	printf("Nverify() AVG took %f milli-seconds \n", (verify_avg / ITERATIONS));

	printf("\nWith M4R in keygen: \n");

	for (int i = 0; i < ITERATIONS; i++) {
		wave_pk_t pk;
		wave_sk_t sk;
		uint8_t salt[SALT_SIZE] = { 0 };

		f3_vector signature = f3_vector_new(N);

		f3_vector m_hash = f3_vector_new(N - K);
		uint8_t mi[MAXMLEN] = { 0 };

		start = clock();
		keygen_m4r(&sk, &pk);
		end = clock();
		timekey_f[i] = (double) (end - start) * 1000.0 / CLOCKS_PER_SEC;
	//	save_matrix(pk, "pk1.bin");

	//	read_matrix(pk, "pk1.bin");
		f3_vector_zero(&signature);
		f3_vector_zero(&m_hash);
		memset(mi, 0, MAXMLEN * sizeof(uint8_t));
		randombytes(mi, MAXMLEN);

		start = clock();
		sign(&signature, &m_hash, salt, mi, MAXMLEN, &sk);
		end = clock();
		time_sign_f[i] = (double) (end - start) * 1000.0 / CLOCKS_PER_SEC;

		size_t alloc = (1 + ((K - 1)) / WORD_LENGTH);
		wave_word a_r0[(1 + ((K - 1)) / WORD_LENGTH)] = { 0 };
		wave_word a_r1[(1 + ((K - 1)) / WORD_LENGTH)] = { 0 };

		memcpy(a_r0, signature.r0 + 45, 88 * sizeof(wave_word));
		memcpy(a_r1, signature.r1 + 45, 88 * sizeof(wave_word));
		shift7(a_r0, alloc);
		shift7(a_r1, alloc);

		start = clock();
		int verify = Nverify(salt, mi, MAXMLEN, a_r0, a_r1, pk);
		end = clock();
		time_verify_f[i] = (double) (end - start) * 1000.0 / CLOCKS_PER_SEC;
		//printf("works? %d\n", verify);

		wave_pk_clear(pk);
		wave_sk_clear(&sk);
		f3_vector_free(&signature);
		f3_vector_free(&m_hash);
	}

	bubbleSort(timekey_f, ITERATIONS - 1);
	bubbleSort(time_verify_f, ITERATIONS - 1);
	bubbleSort(time_sign_f, ITERATIONS - 1);
	sign_avg = 0;
	verify_avg = 0;
	keygen_avg = 0;
	for (int i = 0; i < ITERATIONS; i++) {
		sign_avg += time_sign_f[i];
		verify_avg += time_verify_f[i];
		keygen_avg += timekey_f[i];
	}

	printf("keygen() MEDIAN %f milli-seconds \n", timekey_f[ITERATIONS / 2]);

	printf("keygen() AVG %f milli-seconds \n", (keygen_avg / ITERATIONS));

	printf("sign() MEDIA %f milli-seconds\n", time_sign_f[ITERATIONS / 2]);

	printf("sign() AVG took %f milli-seconds \n", (sign_avg / ITERATIONS));

	printf("Nverify() MEDIAN %f milli-seconds \n",
			time_verify_f[ITERATIONS / 2]);

	printf("Nverify() AVG took %f milli-seconds \n", (verify_avg / ITERATIONS));

	printf("\t FINISHED \n");
	cleanup();
	return EXIT_SUCCESS;
}
