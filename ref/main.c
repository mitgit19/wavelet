#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
