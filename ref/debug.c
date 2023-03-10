#include "debug.h"

void save_matrix(mf3 const *m, char *name){

	FILE *write_ptr;

	write_ptr = fopen(name,"wb");  // w for write, b for binary
	
	for(int i = 0; i < m->n_rows; i++){
		for(int j = 0; j < m->n_cols; j++){
			unsigned char tmp = f3_vector_get_coeff(&m->row[i], j);
			fprintf(write_ptr, &tmp); 
		}
	}
	fclose(write_ptr);

}


void read_matrix(mf3 *m, char* name){

	FILE *read_ptr;
	read_ptr = fopen(name, "rb");
	
	for(int i = 0; i < m->n_rows; i++){
		for(int j = 0; j < m->n_cols; j++){
			unsigned char tmp = 0;
			fscanf(read_ptr, "%c", &tmp);
			f3_vector_setcoeff(&m->row[i], j, tmp);
		}
	}
}

