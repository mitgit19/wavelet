#ifndef COMPRESS_H_
#define COMPRESS_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint-gcc.h>

#include <stdint.h>

#include "vf3.h"
#include "f3.h"

void encode_sig(const f3_vector *sig, const size_t input_len,
		unsigned char *output, size_t *output_len);

void encode_sig_triple(const f3_vector *sig, const size_t input_len,
		unsigned char *output, size_t *output_len);

void encode(const unsigned char *a_h, const unsigned char *a_l, const size_t input_len,
		unsigned char *output, size_t *output_len);

void decode(const unsigned char *input, const size_t input_len,
		unsigned char *a_h,  unsigned char *a_l, size_t *output_len);

void decode_sig(const unsigned char *input, const size_t input_len,
		f3_vector *out);

void decode_sig_triple(const unsigned char *input, const size_t input_len,
		f3_vector *out);


#endif /* COMPRESS_H_ */
