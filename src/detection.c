/* R interface for peak detection function */

#include <R.h>
#include <Rdefines.h>
#include "peak.h"
#include "ncdf.h"

SEXP
peak_detection_main(SEXP Method, SEXP ncData, SEXP Window, SEXP Search,
		    SEXP minInt, SEXP intMatrix)
{
	ncdf_t *nc = !isNull(ncData) ? new_ncdf(ncData) : NULL;
	matrix_t *mat = !isNull(ncData) ? get_intensity_mat(nc) : from_matrix(intMatrix);

	int win     = INTEGER_VALUE(Window);
	int min_int = INTEGER_VALUE(minInt);
	int search  = INTEGER_VALUE(Search);
	const char * method = CHARACTER_VALUE(Method);

	/* create a matrix object, but pass a pointer to pass
	 * the data to be returned to R */
	SEXP ANS = PROTECT(allocMatrix(INTSXP, mat->nr, mat->nc));
	matrix_t *ans = new_mat_alloc(mat->nc, mat->nr, INTEGER_POINTER(ANS));

	peak_detection(method[0], mat, win, search, min_int, ans);

	free_ncdf(nc);
	free_mat(mat);
	free_mat(ans);
	UNPROTECT(1);
	return ANS;
}
