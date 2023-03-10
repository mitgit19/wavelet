#include "hash.h"

/**
 * Function:  expand
 * --------------------
 * expand the input to the correct size in trits.
 * It is algorithm 3 from https://eprint.iacr.org/2021/1432.pdf
 *
 *  input: binary input
 *  len_input: size of input
 *
 * 	h3: structure that will receive the trits
 *
 */

void expand(unsigned char *input, size_t len_input, size_t len_out,
		f3_vector *p) {
	uint8_t buf[168] = { 0 };
	shake256incctx sha3;
	shake256_inc_init(&sha3);

	shake256_inc_absorb(&sha3, input, len_input);
	shake256_inc_squeeze(buf, 168, &sha3);
	size_t len_bytes = 167;
	size_t len = 0;
	while (len < len_out) {
		if (len_bytes == 0) {
			shake256_inc_squeeze(buf, 168, &sha3);
			len_bytes = 167;
		}
		if (buf[len_bytes] < 243) {

			uint8_t to_convert = buf[len_bytes];
			for (int i = 0; i < 8; i = i + 2) {

				uint8_t a = to_convert & 1;
				a = (a | (to_convert & 2)) % 3;
				to_convert /= 3;
				if (len == len_out)
					break;
				f3_vector_setcoeff(p, len, ((a - 1) + 3) % 3);
				//p[len] = ((a - 1) + 3) % 3;
				len++;
			}

		}
		len_bytes--;
	}
	shake256_inc_ctx_release(&sha3);

}

/**
 * Function:  convert_to_trits
 * --------------------
 * convertion from bits to trits
 * It is algorithm 2 from https://eprint.iacr.org/2021/1432.pdf
 *
 *  input: binary input
 *
 * 	h3: structure that will receive the trits
 *
 */

void convert_to_trits(unsigned char *input, f3_vector *h3) {
	size_t index = 0;
	for (int i = 0; i < WORD_LENGTH; i++) {
		uint8_t to_convert = input[i];
		f3_vector_setcoeff(h3, index, (((to_convert % 3) - 1) + 3) % 3);
		to_convert /= 3;
		index++;
	}
}

/*
 * Function:  hash_message
 * --------------------
 * computation of hash message with the vector struct
 *
 *  input:
 *  message: message that we will take the hash
 *  mlen: size of the message
 *
 *  salt: salt that it will be included to hash
 *  msalt: size of the salt
 *
 * 	output:
 *  m_hash: high and low part of the representation of the hash
 *
 *
 *
 */
void hash_message(f3_vector *m_hash, const uint8_t *message, const size_t mlen,
		const uint8_t *salt, const size_t msalt) {
	unsigned char input[MAX_BUFF_SIZE] = { 0 }; // calloc(mlen + msalt, sizeof(unsigned char));
	memcpy(input, salt, msalt * sizeof(uint8_t));
	memcpy(input + msalt, message, mlen * sizeof(uint8_t));

	unsigned char h2[HASH_SIZE] = { 0 };
	sha3_512(h2, input, mlen + msalt);
	f3_vector h3 = f3_vector_new(Mu);
	convert_to_trits(h2, &h3);
	//uint8_t p[N - K - Mu] = { 0 };
	f3_vector p = f3_vector_new(N - K - Mu);
	expand(h2, HASH_SIZE, N - K - Mu, &p);

	memcpy(m_hash->r0, h3.r0, sizeof(wave_word) * h3.alloc);
	memcpy(m_hash->r0 + h3.alloc, p.r0, sizeof(wave_word) * (p.alloc));

	memcpy(m_hash->r1, h3.r1, sizeof(wave_word) * h3.alloc);
	memcpy(m_hash->r1 + h3.alloc, p.r1, sizeof(wave_word) * (p.alloc));

	f3_vector_free(&h3);
	f3_vector_free(&p);

}

/*
 * Function:  hash_message_f
 * --------------------
 * computation of hash message without the vector struct
 *
 *  input:
 *  message: message that we will take the hash
 *  mlen: size of the message
 *
 *  salt: salt that it will be included to hash
 *  msalt: size of the salt
 *
 * 	output:
 *  hash_r0: low part of the representation of the hash
 *  hash_r1: high part of the representation of the hash
 *
 *
 */
void hash_message_f(wave_word *hash_r0, wave_word *hash_r1,
		const uint8_t *message, const size_t mlen, const uint8_t *salt,
		const size_t msalt) {
	unsigned char input[MAX_BUFF_SIZE] = { 0 }; // calloc(mlen + msalt, sizeof(unsigned char));
	memcpy(input, salt, msalt * sizeof(uint8_t));
	memcpy(input + msalt, message, mlen * sizeof(uint8_t));

	unsigned char h2[HASH_SIZE] = { 0 };
	sha3_512(h2, input, mlen + msalt);
	f3_vector h3 = f3_vector_new(Mu);
	convert_to_trits(h2, &h3);
	//uint8_t p[N - K - Mu] = { 0 };
	f3_vector p = f3_vector_new(N - K - Mu);
	expand(h2, HASH_SIZE, N - K - Mu, &p);

	memcpy(hash_r0, h3.r0, sizeof(wave_word) * h3.alloc);
	memcpy(hash_r0 + h3.alloc, p.r0, sizeof(wave_word) * (p.alloc));

	memcpy(hash_r1, h3.r1, sizeof(wave_word) * h3.alloc);
	memcpy(hash_r1 + h3.alloc, p.r1, sizeof(wave_word) * (p.alloc));

	f3_vector_free(&h3);
	f3_vector_free(&p);

}
