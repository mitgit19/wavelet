#ifndef SRC_WAVE_VERIFY_H_
#define SRC_WAVE_VERIFY_H_

#include <stdlib.h>

#include "f3.h"
#include "vf3.h"
#include "mf3.h"

#include "prng.h"

#include "definitions.h"
#include "wave.h"
#include "commons.h"
#include "params.h"
#include "compress.h"

#include "hash.h"

int verify(const uint8_t *salt, const uint8_t *message, const size_t mlen,
		const f3_vector *signature, const mf3 *pk);

int Nverify(const uint8_t *salt, const uint8_t *message, const size_t mlen,
		const wave_word *t_0,const wave_word *t_1 , mf3 *pk);

//int Nverify(f3_vector *hash, const f3_vector *truncated_signature, mf3 *pk);
#endif /* SRC_WAVE_VERIFY_H_ */
