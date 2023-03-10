#ifndef WAVE_H
#define WAVE_H

#include <stdio.h>
#include <stdint.h>

#include "params.h"
#include "definitions.h"




/* Wave parameters, including features for precomputed sampling data */

#define RND_MAX_LEN 16
typedef struct {
    unsigned char size;
    unsigned char value[RND_MAX_LEN];
    prng_t *PRNG;
} random_stream_t;

#define DISTRIB_PREC 4
// represents the floating point number sum(value[i] / 256^(i+1+offset)) >= 0
// if offset >= 0 the number is < 1
typedef struct {
    char offset;
    unsigned char value[DISTRIB_PREC];
} fp_t;

// (cumulated) distribution of size possible values value[i], 0<=i<size
// x[i]-x[i-1] is the probability to get value[i]
// x[size-1] = 1
#define MAX_VALUE_SIZE 569
typedef struct {
    uint16_t size;
    uint16_t value[MAX_VALUE_SIZE];
    fp_t *x;
} Cdist_t;

// array used for rejection
typedef struct {
    int min, max, prec;
    uint8_t **proba;
} farray_t;





#endif
