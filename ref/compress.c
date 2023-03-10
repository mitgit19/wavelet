#include "compress.h"

/*size_t code_len[] = { 2, 2, 3, 3, 4, 4, 4, 5, 5 };

 size_t first_symbol[] = { 0, 0, 0, 2, 4, 7, 10 };

 size_t first_code_r[] = { 0, 0, 0, 4, 12, 30, 32 };

 size_t first_code_l[] = { 0, 0, 0, 16, 24, 30, 32 };

 size_t search_start[] = { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3,
 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5 };*/

size_t code_len[] = { 3, 3, 3, 3, 3, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		7, 9, 10, 10, 10, 10, 10, 10 };

size_t first_symbol[] = { 0, 0, 0, 0, 5, 8, 8, 19, 20, 20, 21, 27 };

size_t first_code_r[] = { 0, 0, 0, 0, 10, 10, 52, 126, 126, 508, 1018, 1024 };

size_t first_code_l[] = { 0, 0, 0, 0, 640, 832, 832, 1008, 1016, 1016, 1018,
		1024 };

size_t search_start[] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 10, 10, 10, 10, 10, 10 };

uint16_t code(size_t pos) {
	return first_code_r[code_len[pos]] + pos - first_symbol[code_len[pos]];
}

uint8_t coeff_isone(const f3_vector *x, int j) {
	uint32_t i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	return (x->r0[i] >> j) & 1;
}

uint8_t coeff_istwo(const f3_vector *x, int j) {
	uint32_t i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	return (x->r1[i] >> j) & 1;
}

void encode_sig(const f3_vector *sig, const size_t input_len,
		unsigned char *output, size_t *output_len) {

	size_t delta = 0;
	size_t byte = 0;
	size_t pos = 0;
	size_t idx = 1;

	for (int i = 0; i < sig->size; i = i + 2) {
		if (coeff_istwo(sig, i)) {
			if (coeff_istwo(sig, i + 1)) {
				pos = 0;
			} else {
				if (coeff_isone(sig, i + 1)) {
					pos = 1;
				} else {
					pos = 4;
				}
			}
		} else {
			if (coeff_istwo(sig, i + 1)) {
				if (coeff_isone(sig, i)) {
					pos = 2;
				} else {
					pos = 6;
				}
			} else {
				if (coeff_isone(sig, i)) {
					if (coeff_isone(sig, i + 1))
						pos = 3;
					else
						pos = 5;
				} else {
					if (coeff_isone(sig, i + 1)) {
						pos = 7;
					} else {
						pos = 8;
					}
				}
			}
		}

		size_t length = code_len[pos];
		delta = delta + length;

		if (delta > 8) {
			unsigned char c = code(pos);
			byte = (byte << length) + code(pos);
			delta = delta - 8;
			output[idx] = (byte >> (delta));
			idx++;
			byte = c & ((1 << delta) - 1);
		} else if (delta == 8) {
			byte = (byte << length) + code(pos);
			output[idx] = byte;
			idx++;
			delta = 0;
			byte = 0;
		} else {
			byte = (byte << length) + code(pos);
		}

	}
	output[0] = delta;
	output[idx] = (byte << (8 - delta));
	idx++;
	*output_len = idx;

}

void encode_sig_triple(const f3_vector *sig, const size_t input_len,
		unsigned char *output, size_t *output_len) {

	size_t delta = 0;
	size_t byte = 0;
	size_t pos = 0;
	size_t idx = 1;

	for (int i = 0; i < sig->size; i = i + 3) {

		if (coeff_istwo(sig, i)) {
			if (coeff_istwo(sig, i + 1)) {
				if (coeff_istwo(sig, i + 2)) {
					pos = 0;
				} else {
					if (coeff_isone(sig, i + 2)) {
						pos = 1;
					} else {
						pos = 8;
					}
				}

			} else {
				if (coeff_istwo(sig, i + 2)) {
					if (coeff_isone(sig, i + 1)) {
						pos = 2;
					} else {
						pos = 9;
					}
				} else {
					if (coeff_isone(sig, i + 1)) {
						if (coeff_isone(sig, i + 2)) {
							pos = 4;
						} else {
							pos = 11;
						}
					} else {
						if (coeff_isone(sig, i + 2)) {
							pos = 12;
						} else {
							pos = 20;
						}
					}

				}
			}

		} else {
			if (coeff_isone(sig, i)) {
				if (coeff_istwo(sig, i + 1)) {
					if (coeff_istwo(sig, i + 2)) {
						pos = 3;
					} else {
						if (coeff_isone(sig, i + 2)) {
							pos = 5;
						} else {
							pos = 14;
						}
					}
				} else {
					if (coeff_istwo(sig, i + 2)) {
						if (coeff_isone(sig, i + 1)) {
							pos = 6;
						} else {
							pos = 15;
						}
					} else {
						if (coeff_isone(sig, i + 1)) {
							if (coeff_isone(sig, i + 2)) {
								pos = 7;
							} else {
								pos = 17;
							}
						} else {
							if (coeff_isone(sig, i + 2)) {
								pos = 18;
							} else {
								pos = 23;
							}
						}
					}
				}
			} else {
				if (coeff_istwo(sig, i + 1)) {
					if (coeff_istwo(sig, i + 2)) {
						pos = 10;
					} else {
						if (coeff_isone(sig, i + 2)) {
							pos = 13;
						} else {
							pos = 21;
						}
					}
				} else {
					if (coeff_istwo(sig, i + 2)) {
						if (coeff_isone(sig, i + 1)) {
							pos = 16;
						} else {
							pos = 22;
						}
					} else {
						if (coeff_isone(sig, i + 1)) {
							if (coeff_isone(sig, i + 2)) {
								pos = 19;
							} else {
								pos = 24;
							}
						} else {
							if (coeff_isone(sig, i + 2)) {
								pos = 25;
							} else {
								pos = 26;
							}
						}

					}
				}
			}

		}

		size_t length = code_len[pos];
		delta = delta + length;
		uint16_t c = code(pos);
		if (delta > 8) {
			if (delta > 16) {
				output[idx] = (byte << (8 - delta + length))
						| (c >> (delta - 8));
				idx++;
				output[idx] = ((c & ((1 << (delta - 8)) - 1)) >> (delta - 16));
				idx++;
				byte = c & ((1 << (delta - 16)) - 1);
				delta = delta - 16;

			} else {
				if (delta == 16) {
					output[idx] = ((byte << (8 - delta + length))
							| (c >> (delta - 8)));
					idx++;
					output[idx] = (c & ((1 << 8) - 1));
					idx++;
					delta = 0;
					byte = 0;

				} else {
					byte = (byte << length) + c;
					output[idx] = (byte >> (delta - 8));
					idx++;
					delta = delta - 8;
					byte = c & ((1 << delta) - 1);

				}

			}

		} else {
			if (delta == 8) {
				byte = (byte << length) + c;
				output[idx] = byte;
				idx++;
				delta = 0;
				byte = 0;
			} else {
				byte = (byte << length) + c;
			}
		}

	}
	output[0] = delta;
	output[idx] = (byte << (8 - delta));
	idx++;
	*output_len = idx;

}

void encode(const unsigned char *a_h, const unsigned char *a_l,
		const size_t input_len, unsigned char *output, size_t *output_len) {

	size_t delta = 0;
	size_t byte = 0;
	size_t pos = 0;
	size_t idx = 1;

	for (int i = 0; i < input_len - 1; i = i + 2) {
		if (a_h[i]) {
			if (a_h[i + 1]) {
				pos = 0;
			} else {
				if (a_l[i + 1]) {
					pos = 1;
				} else {
					pos = 4;
				}
			}
		} else {
			if (a_h[i + 1]) {
				if (a_l[i]) {
					pos = 2;
				} else {
					pos = 6;
				}
			} else {
				if (a_l[i]) {
					if (a_l[i + 1])
						pos = 3;
					else
						pos = 5;
				} else {
					if (a_l[i + 1]) {
						pos = 7;
					} else {
						pos = 8;
					}
				}
			}
		}

		size_t length = code_len[pos];
		delta = delta + length;

		if (delta > 8) {
			unsigned char c = code(pos);
			byte = (byte << length) + code(pos);
			delta = delta - 8;
			output[idx] = (byte >> (delta));
			idx++;
			byte = c & ((1 << delta) - 1);
		} else if (delta == 8) {
			byte = (byte << length) + code(pos);
			output[idx] = byte;
			idx++;
			delta = 0;
			byte = 0;
		} else {
			byte = (byte << length) + code(pos);
		}

	}
	output[0] = delta;
	output[idx] = (byte << (8 - delta));
	idx++;
	*output_len = idx;

}

size_t find_index(size_t buffer) {
	return search_start[buffer];
}

void decode_one(const unsigned char *input, const size_t input_len,
		size_t *index, size_t *delta, size_t *buffer, unsigned char last,
		size_t *s) {
	size_t L = 5;
	size_t l = find_index(*buffer);
	size_t offset = (*buffer - first_code_l[l]) >> (L - l);
	*s = first_symbol[l] + offset;
	size_t length = code_len[*s];

	*buffer = (*buffer & ((1 << (L - length)) - 1)) << length;
	*delta = *delta + length;
	if (*index < input_len) {
		if (*delta <= 8) {
			*buffer += (input[*index] & ((1 << (8 - *delta + length)) - 1))
					>> (8 - *delta);
			return;
		} else {
			if (*index == input_len - 1 && (*delta - 8 - L >= last)) {
				*index = -1;
				*delta = -1;
				*buffer = -1;
				return;

			} else {

				*buffer += (input[*index] & ((1 << (8 - *delta + length)) - 1))
						<< (*delta - 8);
				*buffer += input[*index + 1] >> (16 - *delta);
				*index += 1;
				*delta = *delta - 8;
				return;
			}
		}

	} else {
		if (*index == input_len) {
			if (*delta <= 8) {
				if (*delta < last + L) {
					*buffer += (input[*index]
							& ((1 << (8 - *delta + length)) - 1))
							>> (8 - *delta);
					return;
				} else {
					*index = -1;
					*delta = -1;
					*buffer = -1;
					return;

				}
			} else {
				if (*delta < last + L) {
					*buffer += (input[*index]
							& ((1 << (8 - *delta + length)) - 1))
							<< (*delta - 8);
					*index = *index + 1;
					return;
				} else {
					*index = -1;
					*delta = -1;
					*buffer = -1;
					return;
				}
			}

		} else {
			if (*index == input_len) {
				if (*delta + 8 - L >= last) {
					*index = -1;
					*delta = -1;
					*buffer = -1;
					return;
				} else {
					return;
				}
			}

		}
	}

}

void setvalue_r0(f3_vector *x, size_t j) {
	size_t i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	wave_word z = ((wave_word) 1) << j;
	x->r0[i] |= z;
}
void setvalue_r1(f3_vector *x, size_t j) {
	size_t i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	wave_word z = ((wave_word) 1) << j;
	x->r1[i] |= z;
}

void decode_sig(const unsigned char *input, const size_t input_len,
		f3_vector *out) {
	size_t index = 1;
	size_t L = 5;
	size_t delta = 5;
	size_t buffer = input[1] >> (8 - L);
	unsigned char last = input[0];
	size_t idx = 0;
	size_t s = 0;
	while (index != -1) {

		decode_one(input, input_len, &index, &delta, &buffer, last, &s);
		switch (s) {
		case 0:
			setvalue_r1(out, idx);
			/*	a_h[idx] = 1;
			 a_l[idx] = 0;*/
			idx++;
			setvalue_r1(out, idx);
			/*a_h[idx] = 1;
			 a_l[idx] = 0;*/
			idx++;

			break;
		case 1:
			setvalue_r1(out, idx);
			/*a_h[idx] = 1;
			 a_l[idx] = 0;*/
			idx++;
			/*a_h[idx] = 0;
			 a_l[idx] = 1;*/
			setvalue_r0(out, idx);
			idx++;

			break;

		case 2:
			setvalue_r0(out, idx);
			/*a_h[idx] = 0;
			 a_l[idx] = 1;*/
			idx++;
			/*a_h[idx] = 1;
			 a_l[idx] = 0;*/
			setvalue_r1(out, idx);
			idx++;

			break;

		case 3:
			/*a_h[idx] = 0;
			 a_l[idx] = 1;*/
			setvalue_r0(out, idx);
			idx++;
			/*a_h[idx] = 0;
			 a_l[idx] = 1;*/
			setvalue_r0(out, idx);
			idx++;

			break;

		case 4:
			/*a_h[idx] = 1;
			 a_l[idx] = 0;*/
			setvalue_r1(out, idx);
			idx++;
			/*a_h[idx] = 0;
			 a_l[idx] = 0;*/
			idx++;

			break;

		case 5:
			/*a_h[idx] = 0;
			 a_l[idx] = 1;*/
			setvalue_r0(out, idx);
			idx++;
			/*a_h[idx] = 0;
			 a_l[idx] = 0;*/
			idx++;

			break;

		case 6:
			/*a_h[idx] = 0;
			 a_l[idx] = 0;*/
			idx++;

			/*	a_h[idx] = 1;
			 a_l[idx] = 0;*/
			setvalue_r1(out, idx);
			idx++;

			break;
		case 7:
			/*a_h[idx] = 0;
			 a_l[idx] = 0;*/
			idx++;
			/*a_h[idx] = 0;
			 a_l[idx] = 1;*/
			setvalue_r0(out, idx);
			idx++;

			break;
		case 8:
			/*a_h[idx] = 0;
			 a_l[idx] = 0;*/
			idx++;
			/*a_h[idx] = 0;
			 a_l[idx] = 0;*/
			idx++;
			break;
		default:
			printf("SOmething is wrong\n");
			break;
		}

	}

	out->r0[out->alloc - 1] &= 0x3FFFFFFFFF;
	out->r1[out->alloc - 1] &= 0x3FFFFFFFFF;

}

void decode(const unsigned char *input, const size_t input_len,
		unsigned char *a_h, unsigned char *a_l, size_t *output_len) {

	size_t index = 1;
	size_t L = 5;
	size_t delta = 5;
	size_t buffer = input[1] >> (8 - L);
	unsigned char last = input[0];
	size_t idx = 0;
	size_t s = 0;

	while (index != -1) {

		decode_one(input, input_len, &index, &delta, &buffer, last, &s);
		switch (s) {
		case 0:
			a_h[idx] = 1;
			a_l[idx] = 0;
			idx++;
			a_h[idx] = 1;
			a_l[idx] = 0;
			idx++;

			break;
		case 1:
			a_h[idx] = 1;
			a_l[idx] = 0;
			idx++;
			a_h[idx] = 0;
			a_l[idx] = 1;
			idx++;

			break;

		case 2:
			a_h[idx] = 0;
			a_l[idx] = 1;
			idx++;
			a_h[idx] = 1;
			a_l[idx] = 0;
			idx++;

			break;

		case 3:
			a_h[idx] = 0;
			a_l[idx] = 1;
			idx++;
			a_h[idx] = 0;
			a_l[idx] = 1;
			idx++;

			break;

		case 4:
			a_h[idx] = 1;
			a_l[idx] = 0;
			idx++;
			a_h[idx] = 0;
			a_l[idx] = 0;
			idx++;

			break;

		case 5:
			a_h[idx] = 0;
			a_l[idx] = 1;
			idx++;
			a_h[idx] = 0;
			a_l[idx] = 0;
			idx++;

			break;

		case 6:
			a_h[idx] = 0;
			a_l[idx] = 0;
			idx++;
			a_h[idx] = 1;
			a_l[idx] = 0;
			idx++;

			break;
		case 7:
			a_h[idx] = 0;
			a_l[idx] = 0;
			idx++;
			a_h[idx] = 0;
			a_l[idx] = 1;
			idx++;

			break;
		case 8:
			a_h[idx] = 0;
			a_l[idx] = 0;
			idx++;
			a_h[idx] = 0;
			a_l[idx] = 0;
			idx++;
			break;
		default:
			printf("Something wrong is not correct.\n");
			break;
		}

	}
	*output_len = idx;

}

void decode_one_triple(const unsigned char *input, const size_t input_len,
		size_t *index, size_t *delta, size_t *buffer, unsigned char last,
		size_t *s) {
	size_t L = 10;
	size_t l = find_index(*buffer);
	size_t offset = (*buffer - first_code_l[l]) >> (L - l);
	*s = first_symbol[l] + offset;
	size_t length = code_len[*s];

	*buffer = (*buffer & ((1 << (L - length)) - 1)) << length;
	*delta = *delta + length;

	if (*delta <= 16) {
		if (*index < input_len) {
			if (*delta <= 8) {
				*buffer += (input[*index] & ((1 << (8 - *delta + length)) - 1))
						>> (8 - *delta);
				return;
			} else {
				if (*index == input_len - 1 && *delta - 8 - L >= last) {
					*index = -1;
					*delta = -1;
					*buffer = -1;
					return;
				} else {
					*buffer += (input[*index]
							& ((1 << (8 - *delta + length)) - 1))
							<< (*delta - 8);
					*buffer += input[*index + 1] >> (16 - *delta);
					*index += 1;
					*delta = *delta - 8;
					return;
				}
			}
		} else {
			if (*index == input_len) {
				if (*delta <= 8) {
					if (*delta < last + L) {
						*buffer += (input[*index]
								& ((1 << (8 - *delta + length)) - 1))
								>> (8 - *delta);
						return;
					}
					*index = -1;
					*delta = -1;
					*buffer = -1;
					return;
				} else {
					if (*delta < last + L) {
						*buffer += (input[*index]
								& ((1 << (8 - *delta + length)) - 1))
								>> (8 - *delta);
						*index += 1;
						*delta = *delta - 8;
						return;

					}
					*index = -1;
					*delta = -1;
					*buffer = -1;
					return;

				}

			} else {
				if ((*delta - L + (*index - input_len + 1) * 8) >= last) {
					*index = -1;
					*delta = -1;
					*buffer = -1;
					return;
				}
				return;
			}
		}

	} else { //delta > 16
		if (*index < input_len - 2) {
			*buffer += (input[*index] & ((1 << (8 - *delta + length)) - 1))
					<< (*delta - 8);
			*buffer += input[*index + 1] << (*delta - 16);
			*buffer += input[*index + 2] >> (8 - (*delta - 16));
			*index += 2;
			*delta = *delta - 16;
			return;

		}
		if (*index == input_len - 2) {
			if ((*delta - L - 8) >= last) {
				*index = -1;
				*delta = -1;
				*buffer = -1;
				return;
			}
			*buffer += (input[*index] & ((1 << (8 - *delta + length)) - 1))
					<< (*delta - 8);
			*index += 2;
			*delta = *delta - 16;
			return;

		} else {
			if ((*delta - L - 8) >= last) {
				*index = -1;
				*delta = -1;
				*buffer = -1;
				return;
			}
			*index += 2;
			*delta = *delta - 16;
			return;
		}

	}

}

void decode_sig_triple(const unsigned char *input, const size_t input_len,
		f3_vector *out) {
	size_t index = 2;
	size_t L = 10;
	size_t delta = 2;
	unsigned char last = input[0];
	size_t buffer = (input[1] << (L - 8)) + (input[2] >> (16 - L));

	size_t idx = 0;
	size_t s = 0;
	while (index != -1) {
		decode_one_triple(input, input_len, &index, &delta, &buffer, last, &s);
		switch (s) {

		case 0:
			/*a_h.append(1)
			 a_h.append(1)
			 a_h.append(1)
			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(1)
			 */
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;

			break;
		case 1:
		/*a_h.append(1)
		 a_h.append(1)
		 a_h.append(0)

		 a_l.append(1)
		 a_l.append(1)
		 a_l.append(1)*/
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;
			break;

		case 2:
			/*a_h.append(1)
			 a_h.append(0)
			 a_h.append(1)

			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(1)*/
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;

			break;

		case 3:
			/*a_h.append(0)
			 a_h.append(1)
			 a_h.append(1)
			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(1)*/
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;

			break;

		case 4:
			/*    a_h.append(1)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(1)*/

			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;

			break;

		case 5:
			/*a_h.append(0)
			 a_h.append(1)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(1)*/
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;

			break;

		case 6:
			/*a_h.append(0)
			 a_h.append(0)
			 a_h.append(1)
			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(1)*/
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;

			break;
		case 7:
			/*a_h.append(0)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(1)*/
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;

			break;
		case 8:
			/*a_h.append(1)
			 a_h.append(1)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(0)*/
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;

			idx++;
			break;

		case 9:
			/*a_h.append(1)
			 a_h.append(0)
			 a_h.append(1)
			 a_l.append(1)
			 a_l.append(0)
			 a_l.append(1)*/
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;

			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			break;
		case 10:
			/*a_h.append(0)
			 a_h.append(1)
			 a_h.append(1)
			 a_l.append(0)
			 a_l.append(1)
			 a_l.append(1)*/
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			break;
		case 11:
			/*a_h.append(1)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(0)*/
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;
			idx++;
			break;
		case 12:
			/*a_h.append(1)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(0)
			 a_l.append(1)*/
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			idx++;
			setvalue_r0(out, idx);
			idx++;
			break;
		case 13:
			/*a_h.append(0)
			 a_h.append(1)
			 a_h.append(0)
			 a_l.append(0)
			 a_l.append(1)
			 a_l.append(1)*/
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;
			break;
		case 14:
			/*a_h.append(0)
			 a_h.append(1)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(0)*/
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			idx++;
			break;
		case 15:
			/*a_h.append(0)
			 a_h.append(0)
			 a_h.append(1)
			 a_l.append(1)
			 a_l.append(0)
			 a_l.append(1)*/
			setvalue_r0(out, idx);
			idx++;
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			break;
		case 16:
			/*a_h.append(0)
			 a_h.append(0)
			 a_h.append(1)
			 a_l.append(0)
			 a_l.append(1)
			 a_l.append(1)*/

			idx++;
			setvalue_r0(out, idx);
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			break;
		case 17:
			/*            a_h.append(0)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(1)
			 a_l.append(0)*/
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;
			idx++;
			break;
		case 18:
			/*a_h.append(0)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(0)
			 a_l.append(1)*/
			setvalue_r0(out, idx);
			idx++;
			idx++;
			setvalue_r0(out, idx);
			idx++;
			break;
		case 19:
			/*a_h.append(0)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(0)
			 a_l.append(1)
			 a_l.append(1)*/
			idx++;
			setvalue_r0(out, idx);
			idx++;
			setvalue_r0(out, idx);
			idx++;
			break;
		case 20:
			/*a_h.append(1)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(0)
			 a_l.append(0)*/
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			idx++;
			idx++;
			break;
		case 21:
			/*    a_h.append(0)
			 a_h.append(1)
			 a_h.append(0)
			 a_l.append(0)
			 a_l.append(1)
			 a_l.append(0)*/
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			idx++;
			break;
		case 22:
			/*a_h.append(0)
			 a_h.append(0)
			 a_h.append(1)
			 a_l.append(0)
			 a_l.append(0)
			 a_l.append(1)*/
			idx++;
			idx++;
			setvalue_r1(out, idx);
			setvalue_r0(out, idx);
			idx++;
			break;
		case 23:
			/*a_h.append(0)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(1)
			 a_l.append(0)
			 a_l.append(0)*/
			setvalue_r0(out, idx);
			idx++;
			idx++;
			idx++;
			break;
		case 24:
			/*    a_h.append(0)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(0)
			 a_l.append(1)
			 a_l.append(0)*/
			idx++;
			setvalue_r0(out, idx);
			idx++;
			idx++;
			break;
		case 25:
			/*    a_h.append(0)
			 a_h.append(0)
			 a_h.append(0)
			 a_l.append(0)
			 a_l.append(0)
			 a_l.append(1)*/
			idx++;
			idx++;
			setvalue_r0(out, idx);
			idx++;
			break;
		case 26:
			idx++;
			idx++;
			idx++;
			break;

		default:
			printf("SOmething is wrong\n");
			break;
		}

	}

	out->r0[out->alloc - 1] &= 0x3FFFFFFFFF;
	out->r1[out->alloc - 1] &= 0x3FFFFFFFFF;

}
