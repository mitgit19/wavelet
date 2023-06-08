#include <time.h>
#include <math.h>

#include "cpucycles.h"
#include "vf3.h"
#include "mf3.h"
#include "prng.h"
#include <stdbool.h>

// (WIP) - Temporary section (TODO: Create separate .c and corresponding .h files)

// Deducing the optimal number for k
// size_t m4r_opt_k(size_t n) {
//   size_t k = log2(n);
//   return k;
// }

bool mf3_are_equal(const mf3 *a, const mf3 *b);

// Power 2^k (Considering that for now the implementation is for F2)
size_t two_pow(size_t k);

// MakeTable implementation
void make_table(mf3 *A, size_t r_str, size_t c_str, size_t k, mf3 *T);

// AddRowsFromTable implementation
void add_rows_from_table(mf3 *A, size_t r_str, size_t r_end, size_t c_str, size_t k, mf3 *T);

/* Adding M4RI implementation
Input: A - mxn matrix.
Input: k - an integer k>0. For our case this will be deduced automatically.
Result: A in reduced row echelon form
*/

int m4r_impl(mf3 *A, unsigned int *support, size_t k);