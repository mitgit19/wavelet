#include "m4r.h"

// (WIP) - Temporary section (TODO: Create separate .c and corresponding .h files)

// Deducing the optimal number for k
// size_t m4r_opt_k(size_t n) {
//   size_t k = log2(n);
//   return k;
// }

bool mf3_are_equal(const mf3 *a, const mf3 *b) {

	if ((a->n_cols != b->n_cols) || (a->n_rows != b->n_rows)) {
		return false;
	}

	for (size_t k = 0; k < a->n_rows; k++) {

		if ((a->row[k].size != b->row[k].size) || (a->row[k].alloc != b->row[k].alloc)) {
			return false;
		}
		else {
			for (size_t i = 0; i < a->row[k].alloc; i++) {
				if  (a->row[k].r0[i] != b->row[k].r0[i] || a->row[k].r1[i] != b->row[k].r1[i]) {
					return false;
				}
			}
		}
	}

	return true;
}

// Power 2^k (Considering that for now the implementation is for F2)
size_t two_pow(size_t k) {
  size_t res = 2 << (k-1);
  return res;
}

// MakeTable implementation
void make_table(mf3 *A, size_t r_str, size_t c_str, size_t k, mf3 *T) {

	// TODO: Replace the hardcoded section with separate function
	// Using the pre-processed "fundamental indexes" (the indexes correspond to T | start from 0) 
	size_t I[12];
	memset(I, 0, 12 * sizeof(size_t));

	I[0] = 1;
	I[1] = 2;
	I[2] = 4;
	I[3] = 8;
	I[4] = 16;
	I[5] = 32;
	I[6] = 64;
	I[7] = 128;
	I[8] = 256;
	I[9] = 512;
	I[10] = 1024;
	I[11] = 2048;

	// Setting the first two rows of T after the 0th one
	T->row[1] = A->row[r_str + (k - 1)];

	// Added to ensure that we can work with 1 pivot
	if (k > 1) {
		T->row[2] = A->row[r_str + (k - 2)];
	}

	// TODO: Replace the hardcoded section with separate function
	// Indicator used for keeping track with the upcoming fundamental index
	size_t i = 2;

	// Storing the last fundamental row (to avoid multiple indexing)
	f3_vector fund_row = T->row[2];

	// Storing the number of rows in the prior to the fundamental that will be useful ([0..(I[i]-1)] rule)
	size_t dep_i = 1; // For I[i] = 2 that will be 1

	// Creating T based on 2^k additions (considering min k = 1)
	for (size_t j = 3; j < two_pow(k); j++) {

		// If we are on a fundamental index -> increment by one the indicator
		if (j == I[i]) {
			T->row[j] = A->row[r_str + (k - i - 1)];
			fund_row = T->row[j];
			dep_i = j-1;

			// Ensuring that we won't go out of bounds
			if ((k - 1) > i) {
				i++;
			}
		}

		else {
			// j increments by one but dep_i remains constant thus we effectively start from the row that is dep_row rows above the last fundamental row and move downwards until the row prior to the fundamental itself
			f3_vector_add(&T->row[j - dep_i - 1],&fund_row,&T->row[j]);
		}
	}
}

// AddRowsFromTable implementation
void add_rows_from_table(mf3 *A, size_t r_str, size_t r_end, size_t c_str, size_t k, mf3 *T) {
	
	// Within the set region of rows
	for (size_t r = r_str; r < r_end; r++) {
		size_t index_1 = 0;
		size_t index_2 = 0;

		bool contains_two = false;
		// For the first k columns starting from c_str
		for (size_t j = 0; j < k; j++) {
			uint8_t coeff = f3_vector_get_coeff(&A->row[r], c_str + j);

			// If the coefficient is 2
			if (coeff == 2) {
				contains_two = true;
				coeff = 1;
				index_2 = index_2 + (coeff << (k - j - 1));
			}

			index_1 = index_1 + (coeff << (k - j - 1));
		}

		// Adding the row of T corresponding to the previous index_1 value to the row r of A
		if (contains_two) {
			f3_vector_sub_inplace(&A->row[r], &T->row[index_1]);
			f3_vector_sub_inplace(&A->row[r], &T->row[index_2]);
		}
		else {
			f3_vector_sub_inplace(&A->row[r], &T->row[index_1]);
		}
	}
}

/* Adding M4RI implementation
Input: A - mxn matrix.
Input: k - an integer k>0. For our case this will be deduced automatically.
Result: A in reduced row echelon form
*/

int m4r_impl(mf3 *A, unsigned int *support, size_t k) {

	int pivots[N] = { 0 }; //(int*) malloc(k * sizeof(int));
	int nonpivots[N] = { 0 }; //(int*) malloc(n * sizeof(int));

	int nu_pivs = 0;
	int nu_non_pivs = 0;

	// Initialising both r,c to 0
	size_t r = 0;
	size_t c = 0;

	// Finding number of rows and columns
	size_t m = A->n_rows;
	size_t n = A->n_cols;

	// Finding k based on n_columns (since n>m)
	// size_t k = m4r_opt_k(m);
	// size_t k = 6;

	while((c < n) && (r <= (m - 1))) {

		// If out of boundary (nu of columns)
		if (c+k > n) {
			k = n-c;
		}

		// Computing GaussSubmatrix
		size_t k_bar = m4r_mf3_gauss_elim(A, k, r, c);

		// If pivots found
		if (k_bar > 0) {

			// Storing the pivot positions
			for (size_t j = c; j < c+k_bar; j++) {
				pivots[nu_pivs] = j;
				nu_pivs++;
			}

			// Computing 2^k
			size_t two_pow_k = two_pow(k_bar);

			// Creating table T (size 2^k x n)
			mf3* T = mf3_new(two_pow_k, n);

			// Making the table with all the linear combs for 2^k
			make_table(A,r,c,k_bar,T);

			// TODO: Remove - Used for debugging
			// printf("Creating the table T: \n");
			// mf3_print(T);

			// Adding the rows from the table
			add_rows_from_table(A,0,r,c,k_bar,T);
			add_rows_from_table(A,r+k_bar,m,c,k_bar,T);

		}

		else {
			nonpivots[nu_non_pivs] = c;
			nu_non_pivs++;
		}

		// Updating r,c
		r += k_bar;
		c += k_bar;

		// if (k != k_bar) { // I believe this step is wrongly typed in the M4R algorithm (I saw the GitHub repo has also removed it)
		// 	c += 1;
		// }

		if (k_bar == 0){
			c +=1;
		}
	}

	int pivs_tried = nu_pivs + nu_non_pivs;

	for (size_t j = c; j < n; j++) {
		nonpivots[nu_non_pivs] = j;
		nu_non_pivs++;
	}
	
	memcpy(support, pivots, nu_pivs * sizeof(int));
	memcpy(support + nu_pivs, nonpivots, nu_non_pivs * sizeof(int));

	return pivs_tried;
}