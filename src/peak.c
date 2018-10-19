#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "peak.h"

/** replaces a binary vector of peaks with the respective intensity
 *
 * @param x vector of intensities
 * @param n length of x
 * @param min_int intensiti threshold
 * @param ans a binary vector. if 1, then the value is replaced by the value
 *     of x at the same position if this value is greater than the threshold,
 *     otherwise is set to 0.
 */

void set_max_intensity(int *x, int n, int min_int, int *ans)
{
	for(int j = 0; j < n; j++)
		ans[j] = ans[j] == 1 && x[j] >= min_int ? x[j] : 0;
}

/** performs peak detection based on moving-average smoothing
 *
 * @param x vector of intensities
 * @param win the size of the smoothing windows or the number of points
 *     to average.
 * @param search number of points around a detected peak to find the local
 *     maximum (passed to refine_peak)
 * @param n length of x
 * @param ans a binary vector. if 1, then a peak is found. 0 otherwise
 */
void
smoothing(int *x, int win, int search, int n, int *ans)
{
	double *temp = Calloc(n, double);
	moving(x, win, n, temp);
	int npeaks = find_peak_diff(temp, n, ans);
	refine_peak(x, n, search, ans, npeaks);
	Free(temp);
}

/** performs peak detection based on gaussian smoothing
 *
 * @param x vector of intensities
 * @param coef vector of gaussian coefficients (length = win)
 * @param win the size of the smoothing windows or the number of points
 *     to average.
 * @param search number of points around a detected peak to find the local
 *     maximum (passed to refine_peak)
 * @param n length of x
 * @param ans a binary vector. if 1, then a peak is found. 0 otherwise
 */

void
gaussian (int *x, double *coef, int win, int search, int n, int *ans)
{
	double *temp = Calloc(n, double);
	convolve(x, n, coef, win, temp);
	int npeaks = find_peak_diff(temp, n, ans);
	refine_peak(x, n, search, ans, npeaks);
	Free(temp);
}

/** main function for peak detection
 *
 * Peak detection function. It supports three methods for detecting peaks,
 * ppc, smoothing and gaussian smoothing.
 *
 * @param method a char for selecting the peak picking method. Allowed values
 *     are 'p' (ppc), 's' (smoothing) and ('g) gaussian smoothing.
 * @param mat a pointer to a matrix_t object of intensities. Each column represents a m/z
 *     trace.
 * @param win an integer for the smoothing window ('s' and 'g') or the sliding
 *     window ('p'). Note that in for 'p', the actual window is 2*win+1.
 * @param search an integer for peak search refining. It is only used for
 *     methods 's' and 'g'.
 * @param min_int the minimum intensity threshold
 * @param maxm a pointer to matrix_t which will hold the detected peaks. Same
 *     dimensions of 'mat'
 * @returns an integer. 1 success, 0 fail (invalid method)
 */

int
peak_detection(char method, matrix_t *mat, int win, int search, int min_int,
	       matrix_t *maxm)
{
	double *coef = method == 'g' ? gaussian_coef(win) : NULL;
	int ret = 1;

	for(int k = 0; k < mat->nc; k++) {
		int *x = get_col(mat, k);
		int *ans = get_col(maxm, k);

		switch(method) {
		case 'p':
			peak_detection_ppc(x, win, mat->nr, ans);
			break;
		case 's':
			smoothing(x, win, search, mat->nr, ans);
			break;
		case 'g':
			gaussian(x, coef, win, search, mat->nr, ans);
			break;
		default:
			ret = 0;
			goto clean;
		}

		set_max_intensity(x, mat->nr, min_int, ans);
	}

clean:
	if(coef != NULL)
		Free(coef);
	return ret;
}

