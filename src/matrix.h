#ifndef MATRIX_H
#define MATRIX_H

#include <R.h>

#ifndef USING_R
        #include <stdlib.h>
        #define R_Calloc(n,t) (t *) calloc( n, sizeof(t) )
        #define R_Free(x)     free(x)
	#pragma message "not using R"
#endif

typedef struct
{
	int *x;    /* pointer to data */
	int nc;    /* number of columns */
	int nr;    /* number of rows */
	int mzmin; /* m/z min value */
	int mzmax; /* m/z max value */
	int alloc; /* total memory allocated or 0 */
} matrix_t;

matrix_t * new_mat(int, int);
void free_mat(matrix_t *);

matrix_t * new_mat_empty(void);
matrix_t * new_mat_alloc(int, int, int *);
void mat_add_mz(matrix_t *, int);

#define get_col(mat, k) ( mat->x + ((k) * mat->nr) )

#endif
