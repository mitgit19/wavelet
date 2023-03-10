#ifndef WAVE_YATES_H
#define WAVE_YATES_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "randombytes.h"

#define PARAM_RND_SIZE      16
#define PARAM_RND_BIT_SIZE  128

/**
 *  Generate a 16-bit random number between 0 and `bound-1`
 *
 *  @param[in]  bound  The limit of the number to be generated
 *  @return a 16-bit random number
 **/
uint16_t random_uint16_bounded(uint16_t bound);

/**
 *  Return a uniform random bit
 *
 *  @return random bit 0 or 1
 **/
uint8_t randombit();

/**
 *  Fisher-Yates-Knuth-Yao shuffle on a sequence (in-place)
 *
 *  @note
 *  Fisher-Yates shuffle requires the generation of uniform
 *  random number between a certain range. Knuth-Yao method
 *  may be used to generate such numbers uniformly.
 *
 *  @param[in,out] buffer      The input/output sequence
 **/
void fisher_yates_shuffle(unsigned int *buffer, size_t len);

#endif //WAVE_YATES_H
