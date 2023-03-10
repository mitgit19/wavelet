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

#define ITERATIONS 10

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

int main(void) {

	init(PARAMS_ID);

	printf("Warming up.\n");
	uint8_t compressed[940] = { 0 };
	for (int i = 0; i < ITERATIONS / 2; i++) {
		wave_pk_t pk;
		wave_sk_t sk;
		uint8_t salt[SALT_SIZE] = { 0 };

		memset(compressed, 0, sizeof(uint8_t) * 940);
		f3_vector signature = f3_vector_new(N);

		f3_vector m_hash = f3_vector_new(N - K);
		uint8_t mi[MAXMLEN] = { 0 };

		keygen(&sk, &pk);

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

		f3_vector a = f3_vector_new(K);
		memcpy(a.r0, signature.r0 + 45, 88 * sizeof(wave_word));
		memcpy(a.r1, signature.r1 + 45, 88 * sizeof(wave_word));
		shift7(a.r0, a.alloc);
		shift7(a.r1, a.alloc);
		size_t outlen = 0;

		encode_sig_triple(&a, a.size, compressed, &outlen);

		printf("compressed size: %zu\n", outlen);
		for (int j = 0; j < outlen; j++) {
			printf("%x", compressed[j]);
		}
		printf("\n");
		f3_vector_zero(&a);

		decode_sig_triple(compressed, outlen, &a);

		int res = Nverify(salt, mi, MAXMLEN, a_r0, a_r1, pk);
		if (res != 1) {
			printf("fail...\n");
		}

		wave_pk_clear(pk);
		wave_sk_clear(&sk);
		f3_vector_free(&signature);
		f3_vector_free(&m_hash);
		f3_vector_free(&a);
		printf(".\n");
	}

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
		Nverify(salt, mi, MAXMLEN, a_r0, a_r1, pk);
		end = clock();
		time_verify_f[i] = (double) (end - start) * 1000.0 / CLOCKS_PER_SEC;

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

	cleanup();
	return EXIT_SUCCESS;
}
