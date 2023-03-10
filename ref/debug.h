#ifndef DEBUG_H_
#define DEBUG_H_

#include <stdio.h> 
#include <stdlib.h>

#include "definitions.h"
#include "mf3.h"

void save_matrix(mf3 const *m, char *name);

void read_matrix(mf3 *m, char* name);


#endif
