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
	int *peaks = R_Calloc(np, int);

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
	R_Free(peaks);
}

/* generates gaussian coefficients for gaussian filter (convolution)
 *
 * the coefficients are normalized to abs sum = 1. The resulting vector
 * must be freed in a later called.
 *
 * @param n length of vector. n must be odd.
 */
double * gaussian_coef(int n)
{
	double * coef = R_Calloc(n, double);
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

/**
 * Section for baseline correction
 */

/**
 * returns the k value for quantile computation
 *
 * the relationship between k and p is given by p = (k - 1)/(n - 1)
 * where n is the number of elements. See type 7 in `quantile` man
 * page in R.
 *
 * @param n the number of points of the vector.
 * @param p the probability value in [0, 1] used in the computation of the
 *        quantile.
 * @return the k value such that the formula holds.
 */
static inline int get_k(double n, double p)
{
        return (int) (p * (n - 1) + 1);
}

/**
 * Implements R's quantiles type 7.
 *
 * Computes quantiles in a numeric vector. The vector does not need to be
 * sorted (as previous versions). The sorting is taken care by rPsort.
 *
 * @param xs pointer to a numeric vector
 * @param p the probability 0 <= p <= 1
 * @param n length of the vector
 * @return the computed sample quantile or NAN if error
 * @note NAN's are neither checked for nor handled. It is expected that
 *       there are no NAN's.
 */
double quantile(double *xs, double p, int n)
{
	if(!(p >= 0 && p <= 1))
		return NAN;

	int k = get_k((double) n, p);
	rPsort(xs, n, k - 1);

	/* trivial cases */
	if(p == 0 || p == 1)
		return xs[k - 1];

	rPsort(xs + k, n - k, 0);

	double pk1 = ((double) (k-1)) / ((double) n - 1);
	return xs[k - 1] + (n-1) * (xs[k] - xs[k-1]) * (p - pk1);
}

/**
 * find indices from a to b in a sorted vector.
 *
 * Find indices i* = min(i), j* = max(j), such that
 *    x[i] > a && x[j] < b, if eq == 0
 *    x[i] >= a && x[j] <= b, if eq == 1
 * where i,j = 0, ..., n-1
 *
 * @param x a sorted vector
 * @param a the lower limit
 * @param b the upper limit
 * @param n the length of x
 * @param eq the equivalence parameter (see above)
 * @param pa a pointer to the resulting index pa (see above)
 * @param pb a pointer to the resulting index pb (see above)
 */

void find(double *x, double a, double b, int n, int eq, int *pa, int *pb)
{
	int flag = 0;
	static int ii = 0, jj = 0;

	ii = findInterval(x, n, a, TRUE, FALSE, ii, &flag);
	if(eq == 1 && flag == 0 && a == x[ii - 1])
		ii--;

	jj = findInterval(x, n, b, TRUE, FALSE, jj, &flag);
	if(eq == 0 && x[jj] == b)
		jj--;
	*pa = ii;
	*pb = jj;
}

/**
 * compute sliding quantiles.
 *
 * @param x intensity vector.
 * @param t retention time vector (in arbitrary units, usually seconds).
 * @param qntl the quantile probability to use. 0.5 is the recommended value.
 * @param win the time window around a particular time point (half window to
 *        the left, half window to the right).
 * @param step compute the quantile every `step` steps. This parameter is used
 *        to speed up the computations.
 * @param n the length of vectors `x` and `t`. There is no check that their
 *        length is actually equal.
 * @param ans pointer to where the quantiles will be stored.
 * @return the length of the `ans` vector.
 *
 * @note
 *    It is expected that enough memory has been allocated in the array `ans`.
 *    To be on the safe side, allocate memory to the same length of vectors
 *    `x` and `t`.
 */
int qntl_win(double *x, double *t, double qntl, double win, int step, int n, double *ans)
{
	int a, b, eq = 1, len, k = 0;
	double *tmp = R_Calloc(n, double);

	for(int i = 0; i < n; i+= step) {
		find(t, t[i] - win / 2.0, t[i] + win / 2.0, n, eq, &a, &b);
		len = b - a;
		if(len <= 0) {
			ans[k++] = NAN;
			continue;
		}
		Memcpy(tmp, x + a, len);
		ans[k++] = quantile(tmp, qntl, len);
	}
	R_Free(tmp);
	return k;
}

/**
 * Binary search on a sorted vector
 *
 * @param x pointer to a sorted array (double).
 * @param y the value to search for (double)
 * @param n the length of the array (int).
 * @return (int) the position `i` at which the value `y` can be inserted so that
 *    the order of `x` is preserved. If there is an index `j` such that x[j] == y,
 *    then `j` will be returned.
*/

int binsearch(double *x, double y, int n)
{
	int imin = 0, imax = n, i = n / 2;
	if(n <= 0 || y <= x[0])
		return 0;
	if(y > x[n - 1])
		return n;

	while(imax - imin > 1) {
		if(x[i] == y)
			return i;
		else if(x[i] < y)
			imin = i;
		else
			imax = i;
		i = imin + (imax - imin) / 2;
	}
	return imax;
}
