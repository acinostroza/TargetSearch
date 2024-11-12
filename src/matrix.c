#include "matrix.h"

/**
 * create a new struct matrix with allocated data
 *
 * @param nc number of columns
 * @param nr number of rows
 * @return the allocated matrix struct
 */
matrix_t * new_mat(int nc, int nr)
{
	matrix_t *mat = R_Calloc(1, matrix_t);
	mat->alloc = nc * nr;
	mat->x  = mat->alloc > 0 ? R_Calloc(mat->alloc, int) : NULL;
	mat->nc = nc;
	mat->nr = nr;
	mat->mzmin = 0;
	mat->mzmax = nc - 1;
	return mat;
}

/**
 * create a new empty matrix (no columns or rows)
 *
 * @return the allocated matrix
 */
matrix_t * new_mat_empty(void)
{
	return new_mat(0, 0);
}

/**
 * create a new matrix with pre-allocated data
 *
 * Instead of allocating data, a pointer of pre-allocated data is passed.
 * The number of columns and rows must match the length of the data.
 *
 * @param nc number of columns
 * @param nr number of rows
 * @return the allocated matrix struct
 */
matrix_t * new_mat_alloc(int nc, int nr, int *data)
{
	matrix_t *mat = new_mat_empty();
	mat->nc    = nc;
	mat->nr    = nr;
	mat->mzmin = 0;
	mat->mzmax = nc - 1;
	mat->x     = data;
	mat->alloc = 0;
	return mat;
}

/**
 * sets the mz max and mz min values of a matrix struct
 *
 * @param mz the min m/z.
 */
void mat_add_mz(matrix_t * mat, int mz)
{
	mat->mzmin = mz;
	mat->mzmax = mz + mat->nc + 1;
}

/**
 * free matrix object
 *
 * @param mat the object
 */
void free_mat(matrix_t * mat)
{
	if(mat == NULL)
		return;
	if(mat->alloc != 0)
		R_Free(mat->x);
	R_Free(mat);
}
