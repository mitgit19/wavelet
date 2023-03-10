#ifndef KEYGEN_H_
#define KEYGEN_H_

#include <stdio.h>
#include <stdint.h>

#include "f3.h"
#include "vf3.h"
#include "mf3.h"

#include "prng.h"
#include "randombytes.h"

#include "params.h"
#include "definitions.h"


int keygen(wave_sk_t *sk, wave_pk_t *pk);

void wave_sk_clear(wave_sk_t *sk);

void wave_pk_clear(wave_pk_t pk);

#endif /* KEYGEN_H_ */
