#include "verify.h"
#include <time.h>

int8_t f3_vector_coeff_f(const wave_word *t_0, const wave_word *t_1, int j) {

	int i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	if (((t_0[i] >> j) & 1) && ((t_1[i] >> j) & 1)) {
		return 2;
	} else if (((t_0[i] >> j) & 1)) {
		return 1;
	} else {
		return 0;
	}

}

int Nverify(const uint8_t *salt, const uint8_t *message, const size_t mlen,
		const wave_word *t_0, const wave_word *t_1, mf3 *pk) {

	int cmp = 0;
	size_t alloc = 1 + (((N - K) - 1) / WORD_LENGTH);
	wave_word hash_r0[(1 + (((N - K) - 1) / WORD_LENGTH))] = { 0 };
	wave_word hash_r1[(1 + (((N - K) - 1) / WORD_LENGTH))] = { 0 };
	hash_message_f(hash_r0, hash_r1, message, mlen, salt, SALT_SIZE);

	for (int i = 0; i < K / 2; i++) {
		int8_t c_i, c_ip;

		c_i = (f3_vector_coeff_f(t_0, t_1, 2 * i)
				+ f3_vector_coeff_f(t_0, t_1, (2 * i) + 1)) % 3;
		c_ip = (((f3_vector_coeff_f(t_0, t_1, 2 * i))
				- ((f3_vector_coeff_f(t_0, t_1, 2 * i + 1)))) + 3) % 3;

		if (c_i == 1) {
			f3_vector_sum_inplace_f(hash_r0, hash_r1, pk->row[2 * i].r0,
					pk->row[2 * i].r1, alloc);
		} else if (c_i == 2) {
			f3_vector_sub_inplace_f(hash_r0, hash_r1, pk->row[2 * i].r0,
					pk->row[2 * i].r1, alloc);
		}
		if (c_ip == 1) {
			f3_vector_sum_inplace_f(hash_r0, hash_r1, pk->row[2 * i + 1].r0,
					pk->row[2 * i + 1].r1, alloc);
		} else if (c_ip == 2) {
			f3_vector_sub_inplace_f(hash_r0, hash_r1, pk->row[2 * i + 1].r0,
					pk->row[2 * i + 1].r1, alloc);
		}

	}

	if (f3_vector_coeff_f(t_0, t_1, K - 1) == 1) {
		f3_vector_sum_inplace_f(hash_r0, hash_r1, pk->row[K - 1].r0,
				pk->row[K - 1].r1, alloc);
	} else if (f3_vector_coeff_f(t_0, t_1, K - 1) == 2) {
		f3_vector_sub_inplace_f(hash_r0, hash_r1, pk->row[K - 1].r0,
				pk->row[K - 1].r1, alloc);
	}

	hash_r0[alloc - 1] = hash_r0[alloc - 1] & 0x7F;
	hash_r1[alloc - 1] = hash_r1[alloc - 1] & 0x7F;

	size_t weight_t = popcnt(t_0, 704);
	size_t weight_v = popcnt(hash_r0, alloc * 8);

	if ((weight_t + weight_v) == W)
		cmp = 1;

	//f3_vector_free(&hash);

	return cmp;
}

