#ifndef HASH_H_
#define HASH_H_

#include <stdlib.h>
#include <string.h>

#include "fips202.h"

#include "vf3.h"

#include "prng.h"
#include "randombytes.h"

#include "definitions.h"

#include "commons.h"

#include "wave.h"
#include "params.h"

void hash_message(f3_vector *m_hash, const uint8_t *message, const size_t mlen,
		const uint8_t *salt, const size_t msalt);

void hash_message_f(wave_word *hash_r0, wave_word *hash_r1,
		const uint8_t *message, const size_t mlen, const uint8_t *salt,
		const size_t msalt);

#endif /* HASH_H_ */
