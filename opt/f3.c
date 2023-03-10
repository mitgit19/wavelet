#include "f3.h"
#include "libpopcnt.h"


uint8_t f3_randnonzero(prng_t *PRNG) {
	int8_t c = 1 - (2 * rnd_bit(PRNG)); // 0 or 1
	c = (c + 3) % 3;
	return c; // 1 or -1
}

uint8_t f3_rand(prng_t *PRNG) {
	uint8_t c = rnd_trit(PRNG); // 0, 1, or 2
	return c; // 1, 0, or -1
}

/*
// Hamming weight of v
wave_word f3_array_weight(uint8_t *v, int length) {
	//int i;
	wave_word j = 0;
	j = popcnt(v, length*sizeof(uint8_t));
	return j;
}
*/

// Hamming weight of v
wave_word f3_array_weight(uint8_t *v, int length) {
	//int i;
	wave_word j = 0;
	j = popcnt(v, length*sizeof(uint8_t));
	return j;
}

