
#include "yates.h"


uint8_t randombit()
{
    static int32_t bits_consumed = PARAM_RND_BIT_SIZE;
    static uint8_t rnd_buffer[PARAM_RND_SIZE];
    uint8_t b = 0;

    /**
     * Have we depleted our random source?
     **/
    if (bits_consumed >= PARAM_RND_BIT_SIZE) {
        /**
         * If so, generate PARAM_RND_SIZE bytes
         * of random data as our random source
         **/
        randombytes(rnd_buffer, sizeof(rnd_buffer));
        bits_consumed = 0;
    }

    b = (rnd_buffer[bits_consumed >> 3] & (1 << (bits_consumed & 7))) >> (bits_consumed & 7);
    bits_consumed++;

    return b;
}
uint16_t random_uint16_bounded(uint16_t bound)
{
    uint16_t d, u, x;

    /* Knuth-Yao DDG */
    d = 0; u = 1; x = 0;
    do {
        while (u < bound) {
            u = 2*u;
            x = 2*x + randombit();
        }
        d = u - bound;
        u = d;
    } while (x < d);

    return x - d;
}




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
void fisher_yates_shuffle(unsigned int *buffer, size_t len)
{
    unsigned int index, swap;
    int i = len - 1;
    while (i > 0) {
        index = random_uint16_bounded(i+1);
        swap = buffer[index];
        buffer[index] = buffer[i];
        buffer[i] = swap;
        --i;
    }
}
