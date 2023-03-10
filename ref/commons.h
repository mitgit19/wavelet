#ifndef SRC_UTILS_COMMONS_H_
#define SRC_UTILS_COMMONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "prng.h"

#include "wave.h"

/* identifier, block length, error weight, dimension of U, dimension of V, gap */
/*
 wave_params_t defined_params[] = (wave_params_t []) {
 {"128g", 8492, 7980, 3558, 2047, 81},
 {"96g", 6368, 5984, 2668, 1535, 61},
 {"80g", 5308, 4988, 2224, 1280, 50},
 {"64g", 4246, 3990, 1779, 1024, 40},
 {"", 0, 0, 0, 0, 0} // add new parameters sets above and keep this sentinel row last
 };
 */

/* GLOBAL VARIABLES: precomputed data */
extern farray_t rV, *rU;
extern Cdist_t CLPV, *CLPU;

/* GLOBAL VARIABLES: statistics, others, ... */
extern unsigned long rejU, rejV;

int pickV(prng_t *PRNG);

int acceptV(int t, prng_t *PRNG);

int pickU(int t, prng_t *PRNG);

int acceptU(int t, int l, prng_t *PRNG);

void init(const char *params_id);

void cleanup();

void shift7(wave_word *r, size_t size);

#endif /* SRC_UTILS_COMMONS_H_ */
