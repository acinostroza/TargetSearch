/* peak utilities */

#include "utils.h"

#define _A(x,k,n) ((k) >= 0 && (k) < (n) ? ((double) x[k]) : 0.0)

/**
 * computes the moving average
 * @param x the input array
 * @param L the number of points to average, must be odd
 * @param n the length of the array
 * @param ans the output array. Must be previously allocated.
 */
void moving(const int *x, int L, int n, double *ans)
{
	int w = (L - 1) / 2;
	double K = (double) L;
	/* compute X0 */
	ans[0] = 0;
	for(int i = 0; i <= w; i++)
		ans[0] += _A(x, i, n) / K;

	/* compute the rest recursively */
	for(int i = 0; i < n-1; i++)
		ans[i+1] = ans[i] + (_A(x, i + w + 1, n) - _A(x, i - w, n)) / K;
}

/* convolution of two signals (for smoothing)
 *
 * performs the correlation of two vectors x, y. The output vector is
 * defined as: ans[j] = \sum_{k} x[j - k + o] * y[k], where o is half
 * the length of y. It is required that the length of y is odd.
 *
 * @param x the input array
 * @param n the length of the array
 * @param y the second array
 * @param m length of y. Must be odd.
 * @param ans the result vector. must be previously allocated.
 */
void convolve(const int *x, int n, const double *y, int m, double *ans)
{
	int w = (m - 1) / 2;
	Memzero(ans, n);
	for(int i = 0; i < n; i++)
		for(int k = 0; k < m; k++)
			ans[i] += _A(x, i - w + k, n) * y[k];
}

/* find peaks by change of sign of derivative
 *
 * The output is set to 1 if a change of sign of the first derivative, from
 * positive to negative is detected. set to zero otherwise. The output vector
 * is set to zero and must be previously allocated.
 *
 * @param x the input array
 * @param n the length of x
 * @param ans the resulting array. must be previously allocated.
 * @return number of detected peaks
 */
int find_peak_diff(const double *x, int n, int *ans)
{
	int count = 0;
	Memzero(ans, n);
	for(int i = 1; i < n - 1; i++) {
		if(x[i] > x[i - 1] && x[i] > x[i+1]) {
			ans[i] = 1;
			count++;
			continue;
		}
		/* checks for the unusual case in which there are two maxima
		 * with the same value */
		if(x[i] == x[i+1] && i < n - 2) {
			if(x[i] > x[i-1] && x[i+1] > x[i+2]) {
				ans[i] = 1;
				count++;
			}
		}
	}
	return count;
}

/* search true local maximum around the found (smoothed) maximum
 *
 * The problem with smoothing is that the maximum doesn't
 * always coincide with the true maximum. This function looks for
 * the true maximum around a small search window which should be
 * smaller than the smoothing window.
 *
 * @param x the input array
 * @param n the length of x
 * @param w the search window (left and right)
 * @param np the number of detected peaks (=sum of ans)
 * @param ans a binary of found local maxima (1 maximum, 0 otherwise).
 *        The function refines the maxima.
 */
void refine_peak(const int *x, int n, int w, int *ans, int np)
{
	int k = 0;
	int *peaks = Calloc(np, int);

	/* save detected peaks */
	for(int i = 0; i < n; i++) {
		if(ans[i] == 1) {
			peaks[k++] = i;
		}
		ans[i] = 0;
	}

	for(int i = 0; i < np; i++) {
		int best = peaks[i];
		for(k = peaks[i] - w; k <= peaks[i] + w; k++) {
			if(k < 0 || k >= n)
				continue;
			if(x[k] > x[best]) {
				best = k;
			}
		}
		ans[best] = 1;
	}
	Free(peaks);
}

/* generates gaussian coefficients for gaussian filter (convolution)
 *
 * the coefficiens are normalized to abs sum = 1. The resulting vector
 * must be freed in a later called.
 *
 * @param n length of vector. n must be odd.
 */
double * gaussian_coef(int n)
{
	double * coef = Calloc(n, double);
	double sigma  = (double) (n-1) / 6.0;
	double sum    = 0.0;
	double c      = (n-1) / 2.0;
	for(int i = 0; i < n; i++) {
		coef[i] = exp( -(i - c)*(i - c) / (2*sigma*sigma) );
		sum += coef[i];
	}

	for(int i = 0; i < n; i++)
		coef[i] /= sum;
	return coef;
}

/* This function implements PPC algorithm
 *
 * @param x the input intensity vector
 * @param ispan the search window (left-right). len = 2*ispan+1
 * @param n length of x
 * @param ans output vector of length n. 1 if it is a peak, 0 otherwise
 */
int peak_detection_ppc(int *x, int ispan, int n, int *ans)
{
        int i = ispan, j;
	Memzero(ans, n);

        while (i < n-ispan) {
                ans[i] = 1;
                j = i-ispan;
                while( ans[i] == 1 && j <= i+ispan) {
                        if(x[j] > x[i])
                                ans[i] = 0;
                        j++;
                        if(j == i)
                                j++;

                }
                if(ans[i] == 0)
                        i++;
                if(ans[i] == 1)
                        i += ispan-1;
        }

        return 1;
}

