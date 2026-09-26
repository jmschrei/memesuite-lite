# tomtom.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com> 

import time
import math
import numpy
import numba

from numba import njit
from numba import prange
from numpy import uint64
from numpy import int64


@njit(cache=True)
def _binned_median(x, bins, x_min, x_max, counts):
	"""An internal function for calculating medians quickly.

	This method uses a binning-based approximation to quickly calculate medians
	in linear time with low constants. Rather than using a sorting algorithm,
	which is O(n log n) or more sophisticated approaches that are O(n) but with
	bad constants, this approach approximates the median by dividing the range
	of the array into bins, assigning points to bins in one sweep of the data,
	and then scanning over all bins until half of points have been encountered.

	To get a better approximation of the median, sufficient statistics are
	stored that enable returning the average of all points assigned to the
	bin. When there an odd number of points, or an even number and the middle
	two get assigned to the same bin, this should return the exact median.
	"""

	n, n_bins = len(x), len(bins)
	bins[:] = 0

	halfway = 0
	x_max -= x_min
	for i in range(n):
		z = int((x[i] - x_min) / x_max * (n_bins - 1))
		bins[z, 0] += counts[i]
		bins[z, 1] += x[i] * counts[i]
		halfway += counts[i]

	halfway /= 2
	count = 0
	for i in range(n_bins):
		count += bins[i, 0]
		if count >= halfway:
			return bins[i, 1] / bins[i, 0]
			
	return -99999


@njit(cache=True)
def _binned_median_z(x, bins, x_min, x_max, counts, zb, halfway):
	"""`_binned_median` with the bin indices computed in their own pass.

	The index expression is unchanged, so every index is the same; in its own
	loop it vectorizes. `zb` is scratch of at least len(x) integers. `halfway`
	is sum(counts) / 2, which is the same for every query column and so is
	passed in.

	Only the bin where the scan stops has its values summed. The counts are
	integers, so they are histogrammed as int64 in `bins`'s own memory, and
	the scan finds the stopping bin from them. That bin's sum is then taken
	over the elements that fall in it, in ascending order starting from 0.0,
	which is the same sequence of additions `bins[z, 1] += x[i] * counts[i]`
	made. Blocks of 64 elements with no element in the bin are skipped after
	a check that vectorizes. Bitwise equal to `_binned_median`. `bins` must be
	C-contiguous.
	"""

	n, n_bins = len(x), len(bins)
	x_max -= x_min
	for i in range(n):
		zb[i] = int((x[i] - x_min) / x_max * (n_bins - 1))

	cnt = bins.reshape(-1).view(numpy.int64)[:n_bins]
	cnt[:] = 0
	for i in range(n):
		cnt[zb[i]] += counts[i]

	count = 0
	for b in range(n_bins):
		count += cnt[b]
		if count >= halfway:
			s = 0.0
			n_64 = n - n % 64
			for i0 in range(0, n_64, 64):
				hit = 0
				for i in range(i0, i0 + 64):
					hit |= zb[i] == b
				if hit:
					for i in range(i0, i0 + 64):
						if zb[i] == b:
							s += x[i] * counts[i]

			for i in range(n_64, n):
				if zb[i] == b:
					s += x[i] * counts[i]

			return s / cnt[b]

	return -99999


@njit(cache=True)
def _integer_distances_and_histogram(X, Y, gamma, gamma_int, f, medians, 
	median_bins, X_norm, Y_norm, Y_counts, nq_csum, nq, n_bins):
	"""An internal function for integerized scores and the histogram.

	This function is the main workhorse for the TOMTOM algorithm. It contains
	four conceptual steps: (1) calculate the distance matrix between each column 
	in one query and each column across all targets, (2) subtract out the per
	query-column median, (3) integerize the scores into bins, and (4) calculate
	the histogram of these integers. 

	Several speed efficiencies have been built into this, including caching
	minimum and maximum values for each query column for re-use in the median
	calculation, the binned median approximation, and calculating the histogram
	simultaneously with the binned score matrix. 
	"""
	
	# `gamma` is private scratch, so its buffer is used as (n_rows, n_y): each
	# query column's distances are then contiguous for every loop below.
	n_a, n_y = Y.shape[0], Y.shape[-1]
	g = gamma.reshape((gamma.shape[1], gamma.shape[0]))
	x2 = numpy.empty(n_a, dtype=numpy.float64)
	zb = numpy.empty(n_y, dtype=numpy.int32)
	mxb = numpy.empty(64, dtype=numpy.float64)
	mnb = numpy.empty(64, dtype=numpy.float64)
	halfway = 0
	for j in range(n_y):
		halfway += Y_counts[j]
	halfway /= 2

	# Calculate the Euclidean distance between query and targets. The terms are
	# subtracted in the same order as before; `2 * X` is hoisted, which is the
	# same product. The min/max pass is separate so the distance loop
	# vectorizes. It runs as 64 lanes through numpy.maximum/minimum, which
	# compile to packed code where a serial max/min chain cannot. Every g value
	# is -sqrt(z) with z > 0 or +0.0: never NaN and never -0.0, so max and min
	# give the same bits in any order.
	z_min, z_max = 9999999.9, -9999999.9
	for i in range(nq):
		z_min_, z_max_ = 9999999.9, -9999999.9
		xn = X_norm[i + nq_csum]

		if n_a == 4:
			x0 = 2 * X[0, i + nq_csum]
			x1 = 2 * X[1, i + nq_csum]
			x2_ = 2 * X[2, i + nq_csum]
			x3 = 2 * X[3, i + nq_csum]
			for j in range(n_y):
				z = xn + Y_norm[j]
				z -= x0 * Y[0, j]
				z -= x1 * Y[1, j]
				z -= x2_ * Y[2, j]
				z -= x3 * Y[3, j]
				g[i, j] = -math.sqrt(z) if z > 0 else 0
		else:
			for k in range(n_a):
				x2[k] = 2 * X[k, i + nq_csum]
			for j in range(n_y):
				z = xn + Y_norm[j]
				for k in range(n_a):
					z -= x2[k] * Y[k, j]
				g[i, j] = -math.sqrt(z) if z > 0 else 0

		nb = n_y - n_y % 64
		if nb > 0:
			mxb[:] = g[i, :64]
			mnb[:] = g[i, :64]
			for j in range(64, nb, 64):
				numpy.maximum(mxb, g[i, j:j+64], mxb)
				numpy.minimum(mnb, g[i, j:j+64], mnb)
			for u in range(64):
				z_max_ = max(z_max_, mxb[u])
				z_min_ = min(z_min_, mnb[u])
		for j in range(nb, n_y):
			z = g[i, j]
			z_max_ = max(z_max_, z)
			z_min_ = min(z_min_, z)
		
		# Subtract out the median from each row
		m = _binned_median_z(g[i], median_bins, z_min_, z_max_, Y_counts,
			zb, halfway)
		medians[i] = m

		z_min = min(z_min, z_min_ - m)
		z_max = max(z_max, z_max_ - m)
			
	# Find the minimum value and the number of bins needed to get there
	i_min = int(math.floor(z_min)) #offset
	bin_scale = int(math.floor(n_bins / (z_max - i_min))) #scale
	offset = -i_min * bin_scale

	for i in range(nq):
		medians[i] = medians[i] + i_min
	
	f[:] = 0
	ys = numpy.sum(Y_counts)
	w = numpy.empty(n_y, dtype=numpy.float64)
	for j in range(n_y):
		w[j] = Y_counts[j] / ys

	# Convert the distances to bins and record the histogram of counts. The
	# bin indices are computed in their own loop, which vectorizes, and the
	# scatter-add then runs in the original order. Every index lies in
	# [0, n_bins], because f has n_bins + 1 columns, so int32 holds it exactly.
	zb = numpy.empty(n_y, dtype=numpy.int32)
	for i in range(nq):
		k = nq - i - 1
		mi = medians[i]
		for j in range(n_y):
			zb[j] = math.floor((g[i, j] - mi) * bin_scale + 0.5)

		for j in range(n_y):
			x = zb[j]
			gamma_int[j, k] = x - offset
			f[i, uint64(x)] += w[j]

	return uint64(offset)


@njit(cache=True)
def _pairwise_max(x, y, y_csum, z, n):
	"""An internal function for the pdf of the maximum of two pdfs.

	This function takes in two probability distribution functions and
	returns the probability distribution function for the maximum of
	the two. In other words, it returns the probability distribution
	for the maximum of a randomly drawn sample from the first
	distribution and a randomly drawn sample from the second
	distribution.

	This function relies on knowing that the cumsum of y will be
	precalculated and that the cumsum of x has to be recalculated
	each call.
	"""
	
	if x[0] == -1:
		z[:] = y[:]        
	else:
		x_csum = 0
		for i in range(n):
			x_csum += x[i]
			z[i] = x[i] * y_csum[i] + y[i] * x_csum - x[i] * y[i]

 
@njit(cache=True)
def _pairwise_max_window(x, y, y_csum, z, L, H, copy):
	"""`_pairwise_max` restricted to the bins [L, H).

	x and y must be zero outside [L, H), and so z is too; only [L, H) of z is
	written. When `copy` is true, x is the empty starting point and z = y. The
	running sum of x starts at L, where it is still exactly zero.
	"""

	if copy:
		for i in range(L, H):
			z[i] = y[i]
	else:
		x_csum = 0.0
		for i in range(L, H):
			x_csum += x[i]
			z[i] = x[i] * y_csum[i] + y[i] * x_csum - x[i] * y[i]


@njit(cache=True, inline='always')
def _pairwise_max_support(x, y, y_csum, z, L, H, x_lo, x_hi, y_lo, y_hi, 
	copy, inplace):
	"""`_pairwise_max_window` using the nonzero supports of x and y.

	Within [L, H), x is zero outside [x_lo, x_hi) and y outside [y_lo, y_hi),
	and y's support is not empty. Then y_csum is zero below y_lo, the running
	sum of x is zero below x_lo and constant from x_hi, and every term of
	`x[i] * y_csum[i] + y[i] * x_csum - x[i] * y[i]` whose factor is zero is
	exactly +0.0. So z is exactly zero below max(x_lo, y_lo) and above
	max(x_hi, y_hi), equals `x[i] * y_csum[i]` where only x is nonzero and
	`y[i] * x_csum` where only y is, and elsewhere is computed by the full
	expression with the same running sum: skipped additions to it are all of
	+0.0. z is written over all of [L, H), as `_pairwise_max_window` does, 
	except that an in-place z (z is x) is already zero outside [x_lo, x_hi).
	Returns the support of z.
	"""

	if copy:
		for i in range(L, y_lo):
			z[i] = 0.0
		for i in range(y_lo, y_hi):
			z[i] = y[i]
		for i in range(y_hi, H):
			z[i] = 0.0
		return y_lo, y_hi

	z_lo, z_hi = max(x_lo, y_lo), max(x_hi, y_hi)
	m = min(x_hi, y_hi)

	x_csum = 0.0
	for i in range(x_lo, min(z_lo, x_hi)):
		x_csum += x[i]

	if inplace:
		for i in range(x_lo, min(z_lo, x_hi)):
			z[i] = 0.0
	else:
		for i in range(L, z_lo):
			z[i] = 0.0
		for i in range(z_hi, H):
			z[i] = 0.0

	for i in range(z_lo, m):
		x_csum += x[i]
		z[i] = x[i] * y_csum[i] + y[i] * x_csum - x[i] * y[i]

	if x_hi > y_hi:
		for i in range(max(y_hi, z_lo), x_hi):
			z[i] = x[i] * y_csum[i]
	else:
		for i in range(max(x_hi, z_lo), y_hi):
			z[i] = y[i] * x_csum

	return z_lo, z_hi


@njit(cache=True, inline='always')
def _pm(x, A, A_csum, a, b, z, L, H, a_lo, a_hi, x_lo, x_hi, copy, inplace):
	"""One step of the B build: z = max(x, A[a, b]) over [L, H), returning 
	the support of z. x_lo and x_hi are x's support; unused when `copy`."""

	if a_lo[a, b] < 0:
		_pairwise_max_window(x, A[a, b], A_csum[a, b], z, L, H, copy)
		return L, H

	return _pairwise_max_support(x, A[a, b], A_csum[a, b], z, L, H, x_lo, 
		x_hi, a_lo[a, b], a_hi[a, b], copy, inplace)


@njit(cache=True, inline='always')
def _pm_fused2(x, y0, c0, z0, y1, c1, z1, lo, hi):
	"""Two chained steps of `_pairwise_max_window` in one pass over [lo, hi):
	z0 = max(x, y0), then z1 = max(z0, y1).

	Each step is the same expression with its own running sum, accumulated in
	the same order, so each z is bitwise what two separate passes give. The
	two running sums overlap in one pass. x may be z0 or z1 (in place): x[i] 
	is read before z0[i] or z1[i] is written.
	"""

	s0, s1 = 0.0, 0.0
	for i in range(lo, hi):
		v = x[i]
		s0 += v
		t0 = v * c0[i] + y0[i] * s0 - v * y0[i]
		z0[i] = t0
		s1 += t0
		z1[i] = t0 * c1[i] + y1[i] * s1 - t0 * y1[i]


@njit(cache=True, inline='always')
def _pm_chain(B, src, sa, sb, dst, n_steps, A, A_csum, L, H, a_lo, a_hi, 
	x_lo, x_hi):
	"""Run n_steps steps of the B build: step k sets B[dst[k]] to the max of
	the previous result (B[src] for the first) and A[sa[k], sb[k]]. Returns 
	the support of the last result.

	Steps are fused two at a time with `_pm_fused2`, over [lo, hi): lo is 
	the start of the input's support, and hi the largest end of the input's 
	and every A row's support. Every result is supported in [lo, hi) (all of
	[L, H) for an all-zero A row, whose window kernel runs over [L, H)), and 
	outside the supports every term of the expression is exactly +0.0, so 
	the full expression over [lo, hi) gives the bits `_pm` gives. A result 
	row that is not the input row is zeroed over the rest of [L, H), as a 
	non-in-place `_pm` does. The remaining steps go through `_pm`.
	"""

	k, row = 0, src
	while k + 2 <= n_steps:
		lo, hi, z_lo, z_hi = x_lo, x_hi, x_lo, x_hi
		for r in range(k, k+2):
			a, b = sa[r], sb[r]
			if a_lo[a, b] < 0:
				lo, hi, z_lo, z_hi = min(lo, L), H, L, H
			else:
				hi = max(hi, a_hi[a, b])
				z_lo, z_hi = max(z_lo, a_lo[a, b]), max(z_hi, a_hi[a, b])

		_pm_fused2(B[row], A[sa[k], sb[k]], A_csum[sa[k], sb[k]], B[dst[k]],
			A[sa[k+1], sb[k+1]], A_csum[sa[k+1], sb[k+1]], B[dst[k+1]], lo, hi)

		for r in range(k, k+2):
			d = dst[r]
			if d != row and (r == k or d != dst[r-1]):
				for i in range(L, lo):
					B[d, i] = 0.0
				for i in range(hi, H):
					B[d, i] = 0.0

		row, x_lo, x_hi = dst[k+1], z_lo, z_hi
		k += 2

	while k < n_steps:
		x_lo, x_hi = _pm(B[row], A, A_csum, sa[k], sb[k], B[dst[k]], L, H, 
			a_lo, a_hi, x_lo, x_hi, False, dst[k] == row)
		row = dst[k]
		k += 1

	return x_lo, x_hi


@njit(cache=True, inline='always')
def _A_cumsum_fill(A_csum, i, j, lo, hi, n_bins, c, n, acc):
	"""The constant parts of the row A_csum[i, j] around its running sum."""

	A_csum[i, j, :lo] = 0
	A_csum[i, j, hi+1:n_bins*(j+1)+c] = acc
	A_csum[i, j, n_bins*(j+1)+c:n] = 1


@njit(cache=True)
def _A_cumsum(A, A_csum, nq, n_bins, offset, n):
	"""An internal function for the cumulative sums of the span backgrounds.

	A[i, j] can only be nonzero in [c + m, c + m*n_bins], where m = j - i + 1
	is the span length, so the running sum is 0.0 below that range and the
	total above it. Up to n_bins*(j+1) + c it holds the total and past that
	1, over the first n entries, which are all that `_pairwise_max` reads.

	c depends only on m, so every row of one span length has the same range.
	Two such rows are summed in one pass, each in its own accumulator, so
	the two serial chains overlap. Each row's sum is the same chain of 
	additions in the same order as a pass over that row alone.
	"""

	for m in range(1, nq+1):
		m = uint64(m)
		c = uint64(offset * (nq - m))
		lo, hi = c + m, c + m*n_bins

		i = 0
		while i + 2 <= nq - m + 1:
			i0 = uint64(i)
			a0, s0 = A[i0, i0+m-1], A_csum[i0, i0+m-1]
			a1, s1 = A[i0+1, i0+m], A_csum[i0+1, i0+m]

			acc0, acc1 = 0.0, 0.0
			for k in range(lo, hi+1):
				acc0 += a0[k]
				acc1 += a1[k]
				s0[k] = acc0
				s1[k] = acc1

			_A_cumsum_fill(A_csum, i0, i0+m-1, lo, hi, n_bins, c, n, acc0)
			_A_cumsum_fill(A_csum, i0+1, i0+m, lo, hi, n_bins, c, n, acc1)
			i += 2

		while i < nq - m + 1:
			i0 = uint64(i)
			a0, s0 = A[i0, i0+m-1], A_csum[i0, i0+m-1]
			acc0 = 0.0
			for k in range(lo, hi+1):
				acc0 += a0[k]
				s0[k] = acc0

			_A_cumsum_fill(A_csum, i0, i0+m-1, lo, hi, n_bins, c, n, acc0)
			i += 1


@njit(cache=True)
def _convolve_span(prev, row, f, k_lo, k_hi, prev_base, row_base):
	"""An internal function for one step of the span convolution.

	Adds prev[k + prev_base] * f[s] into row[k + row_base + s] for every k in
	[k_lo, k_hi] and s in [0, len(f)). Each element of row receives its terms
	in ascending k, exactly as a scalar scatter over k would add them.

	With a = prev[k_lo + prev_base:] and out = row[k_lo + row_base:], this is
	out[q] = sum over k' of a[k'] * f[q - k']. For a fixed q, ascending k' is
	descending s = q - k', so the passes run over s from high to low, each
	adding f[s] * a[t] into out[t + s]. The vector then runs over a, which is
	several times longer than f. Four s are taken per pass, highest first, so
	each element of out is loaded and stored once per four terms. Each term is
	the same product (f[s] * a == a * f[s]) added in the same order, so the
	result is bitwise the same. A term with a zero factor, which the scatter
	skipped, adds +0.0 to a non-negative value and changes nothing.
	"""

	L = f.shape[0]
	nk = k_hi - k_lo + 1

	if nk < 4 or L < 4:
		for k in range(k_lo, k_hi+1):
			a = prev[k+prev_base]
			if a != 0:
				dst = row[k+row_base:k+row_base+L]
				for s in range(L):
					dst[s] += a * f[s]
		return

	a = prev[k_lo+prev_base:k_hi+prev_base+1]
	out = row[k_lo+row_base:k_hi+row_base+L]
	n = nk - 3
	v0, v1, v2, v3 = a[0:n], a[1:1+n], a[2:2+n], a[3:3+n]

	s0 = L - 1
	while s0 >= 3:
		c0, c1, c2, c3 = f[s0], f[s0-1], f[s0-2], f[s0-3]
		b = s0 - 3

		# The first and last three outputs lack some of the four terms.
		for e in range(3):
			x = out[b+e]
			for r in range(3-e, 4):
				x += f[s0-r] * a[e-3+r]
			out[b+e] = x

		# Contiguous 1-D views keep every index a non-negative loop
		# counter, which is what lets this loop vectorize.
		dst = out[b+3:b+3+n]
		for t in range(n):
			x = dst[t]
			x += c0 * v0[t]
			x += c1 * v1[t]
			x += c2 * v2[t]
			x += c3 * v3[t]
			dst[t] = x

		for e in range(3):
			x = out[b+nk+e]
			for r in range(3-e):
				x += f[s0-r] * a[n+e+r]
			out[b+nk+e] = x

		s0 -= 4

	while s0 >= 0:
		c = f[s0]
		if c != 0:
			dst = out[s0:s0+nk]
			for t in range(nk):
				dst[t] += c * a[t]
		s0 -= 1


@njit(cache=True)
def _p_value_backgrounds(f, A, B, A_csum, nq, n_bins, t_max, offset, 
	needed=None):
	"""An internal function that calculates the backgrounds for p-values.

	This method takes in the histogram of integerized scores `f` and returns 
	the background probabilities of each overlap achieving a given score. 
	These scores are calculated for the complete overlap of the query and
	target, but also for all overhangs where only part of the query and the
	target are overlapping (on either end). Additionally, background
	probabilities are calculated for all spans across the query for when the
	target is smaller than the query and has to be scanned against it.

	When `needed` is given, only the rows `t` with `needed[t]` true are
	finished; the others hold unspecified values.
	"""

	n = n_bins*nq + nq*offset

	# First and last nonzero bin of each query column's histogram. Every term
	# is non-negative, so a skipped zero term would only have added +0.0, and
	# the remaining terms are still added in the same order: bitwise-exact.
	f_lo = numpy.empty(nq, dtype='int64')
	f_hi = numpy.empty(nq, dtype='int64')
	for j in range(nq):
		f_lo[j], f_hi[j] = n_bins+1, 0
		for l in range(1, n_bins+1):
			if f[j, l] != 0:
				f_hi[j] = l
				if f_lo[j] > n_bins:
					f_lo[j] = l
	
	for i in range(nq):
		i = uint64(i)

		# Bounds on the nonzero support of A[i, j-1], in the `k` coordinate
		k_lo, k_hi = 0, -1
		for j in range(i, nq):
			j, c = uint64(j), uint64(offset * (nq - j + i - 1))

			# Only rows A[i, j] with i <= j < nq are used, so clear those
			# rather than the whole Q_max x Q_max x n_len workspace. The whole
			# row is cleared, not just the first n bins that are read, so that
			# each row still holds a complete distribution.
			A[i, j] = 0
			
			if i == j:
				for l in range(1, n_bins+1):
					l = uint64(l)
					A[i, j, l+c] = f[j, l]

				k_lo, k_hi = f_lo[j], f_hi[j]
			else:
				l_lo, l_hi = f_lo[j], f_hi[j]

				_convolve_span(A[i, j-1], A[i, j], f[j, l_lo:l_hi+1], max(k_lo, 0),
					min(k_hi, numpy.int64(n_bins*j)), int64(c+offset), int64(c)+l_lo)

				k_lo, k_hi = k_lo + l_lo, k_hi + l_hi


	_A_cumsum(A, A_csum, nq, n_bins, offset, n)

	###

	# Every A[i, j] is zero outside [lo, hi), where lo and hi follow from the
	# lowest and highest nonzero bins of f, so every B row is too: the pdf of
	# a maximum is zero below the larger of the two lower ends and above the
	# larger of the two upper ends. B is built only over [L, H), the union of
	# those windows. The values dropped are exact zeros, so every value that
	# is kept is computed by the same operations in the same order.
	L, H = int64(n), int64(0)
	a_lo = numpy.empty((nq, nq), dtype='int64')
	a_hi = numpy.empty((nq, nq), dtype='int64')
	for i in range(nq):
		lo, hi, empty = int64(0), int64(0), False
		for j in range(i, nq):
			if f_lo[j] > n_bins:
				empty = True
			if empty:
				a_lo[i, j] = -1
				continue

			lo += f_lo[j]
			hi += f_hi[j]
			c = int64(offset) * int64(nq - j + i - 1)
			a_lo[i, j], a_hi[i, j] = lo + c, hi + c + 1
			L = min(L, lo + c)
			H = max(H, hi + c + 1)

	if H > int64(n):
		H = int64(n)
	if H <= L:
		L, H = 0, 0

	# A[a, b] is zero outside [a_lo, a_hi) (clipped to [L, H)), and a_lo is 
	# -1 when the row is all zero; its A_csum is then not zero below the fill
	# of 1, so it goes through the unrestricted `_pairwise_max_window` and 
	# the result's support is taken to be all of [L, H).
	for i in range(nq):
		for j in range(i, nq):
			if a_lo[i, j] >= 0:
				a_hi[i, j] = min(a_hi[i, j], H)

	# Rows below nq are rebuilt from A alone further down, so the first pass
	# only matters as the start of the chain that builds rows nq..t_max. 
	# Each row is built by a chain of steps that `_pm_chain` runs.
	n_steps_max = 2*nq + t_max
	sa = numpy.empty(n_steps_max, dtype='int64')
	sb = numpy.empty(n_steps_max, dtype='int64')
	dst = numpy.empty(n_steps_max, dtype='int64')

	if H > L and t_max >= nq:
		b_lo, b_hi, ns = L, H, 0
		if nq > 1:
			b_lo, b_hi = _pm(B[0], A, A_csum, 0, 0, B[1], L, H, a_lo, a_hi,
				b_lo, b_hi, True, False)
			sa[ns], sb[ns], dst[ns] = nq-1, nq-1, 1
			ns += 1

		for i in range(2, nq):
			sa[ns], sb[ns], dst[ns] = 0, i-1, i
			sa[ns+1], sb[ns+1], dst[ns+1] = nq-i, nq-1, i
			ns += 2

		for i in range(nq, t_max+1):
			if i == 1:
				b_lo, b_hi = _pm(B[0], A, A_csum, 0, nq-1, B[1], L, H, a_lo,
					a_hi, b_lo, b_hi, True, False)
				continue
			sa[ns], sb[ns], dst[ns] = 0, nq-1, i
			ns += 1

		b_lo, b_hi = _pm_chain(B, 1, sa, sb, dst, ns, A, A_csum, L, H, a_lo,
			a_hi, b_lo, b_hi)

	for i in range(1, min(nq, t_max+1)):
		if H <= L:
			break
		if needed is not None and not needed[i]:
			continue

		b_lo, b_hi = _pm(B[i], A, A_csum, 0, i-1, B[i], L, H, a_lo, a_hi, 
			L, H, True, True)

		ns = 0
		for j in range(1, nq - i + 1):
			sa[ns], sb[ns], dst[ns] = j, j+i-1, i
			ns += 1
	
		for j in range(i-1):
			sa[ns], sb[ns], dst[ns] = 0, j, i
			sa[ns+1], sb[ns+1], dst[ns+1] = nq-1-j, nq-1, i
			ns += 2

		b_lo, b_hi = _pm_chain(B, i, sa, sb, dst, ns, A, A_csum, L, H, a_lo,
			a_hi, b_lo, b_hi)

	# Row 0 is the all -1 starting point and is never a real distribution.
	if needed is None or needed[0]:
		for j in range(n):
			B[0, j] = -1
		for j in range(1, n):
			B[0, j] += B[0, j-1]
		for j in range(n):
			b = 1 - B[0, j]
			B[0, j] = b if b > 0 else 0.0

	# `axis` is not implemented for cumsum. The pdf is zero below L, so the
	# CDF is zero there and the survival 1, and it is flat from H on.
	#
	# The cumsum accumulates floating-point round-off across thousands of bins
	# and the underlying distribution does not sum to exactly 1, so at the
	# extreme right tail (the very best matches) the CDF can land just above
	# 1.0. A survival probability cannot be negative, so clamp the round-off
	# to zero to avoid returning tiny negative p-values.
	rows = numpy.empty(t_max+1, dtype=numpy.int64)
	nr = 0
	for i in range(1, t_max+1):
		if needed is None or needed[i]:
			rows[nr] = i
			nr += 1

	# Each row's running sum is its own serial chain of additions, in the same
	# order as `B[i, j] += B[i, j-1]` (a + b == b + a bitwise). Four rows are
	# summed at once so that four independent chains overlap their latency,
	# and each survival value is computed from the sum as it is produced.
	r = 0
	while r + 4 <= nr:
		b0, b1, b2, b3 = B[rows[r+0]], B[rows[r+1]], B[rows[r+2]], B[rows[r+3]]
		if H > L:
			a0, a1, a2, a3 = b0[L], b1[L], b2[L], b3[L]
			s = 1 - a0
			b0[L] = s if s > 0 else 0.0
			s = 1 - a1
			b1[L] = s if s > 0 else 0.0
			s = 1 - a2
			b2[L] = s if s > 0 else 0.0
			s = 1 - a3
			b3[L] = s if s > 0 else 0.0
			for j in range(L+1, H):
				a0 = b0[j] + a0
				a1 = b1[j] + a1
				a2 = b2[j] + a2
				a3 = b3[j] + a3
				s0 = 1 - a0
				b0[j] = s0 if s0 > 0 else 0.0
				s1 = 1 - a1
				b1[j] = s1 if s1 > 0 else 0.0
				s2 = 1 - a2
				b2[j] = s2 if s2 > 0 else 0.0
				s3 = 1 - a3
				b3[j] = s3 if s3 > 0 else 0.0
		r += 4

	while r < nr:
		b0 = B[rows[r]]
		if H > L:
			a0 = b0[L]
			s = 1 - a0
			b0[L] = s if s > 0 else 0.0
			for j in range(L+1, H):
				a0 = b0[j] + a0
				s = 1 - a0
				b0[j] = s if s > 0 else 0.0
		r += 1

	for r in range(nr):
		i = rows[r]
		for j in range(L):
			B[i, j] = 1.0
		tail = B[i, H-1] if H > 0 else 1.0
		for j in range(H, n):
			B[i, j] = tail
			

@njit(cache=True)
def _p_values(gamma, B_cdfs, rr_inv, T_lens, iq, nq, offset, results,
	reverse_complement=1):
	"""An internal function for calculating the best match and p-values.

	Chooses the width of the running sums and calls `_p_values_sums`. Each
	sum is nq * offset plus at most nq int16 values of `gamma`, so it lies in
	[-nq * 32768, nq * (offset + 32767)]. When that fits in int32 the sums
	are exact in int32, which halves the accumulation's vector width;
	otherwise they are kept in int64.
	"""

	# Sized by the longest target, not by gamma, whose rows are the unique
	# target columns and can be fewer than a target's length after hashing.
	n_sums = T_lens.max() + nq - 1

	if int64(nq) * (int64(offset) + 32768) <= 2147483647:
		t_sums = numpy.empty(n_sums, dtype='int32')
		_p_values_sums(gamma, B_cdfs, rr_inv, T_lens, iq, nq, offset, 
			results, reverse_complement, t_sums)
	else:
		t_sums64 = numpy.empty(n_sums, dtype='int64')
		_p_values_sums(gamma, B_cdfs, rr_inv, T_lens, iq, nq, offset, 
			results, reverse_complement, t_sums64)


@njit(cache=True)
def _p_values_sums(gamma, B_cdfs, rr_inv, T_lens, iq, nq, offset, results,
	reverse_complement, t_sums):
	"""An internal function for calculating the best match and p-values.

	This function will take in the integerized score matrix `gamma` and
	background distributions `B_cdfs` and calculate the best overlap.
	The best overlap is calculated as the best sum of scores across the
	alignment, minus a penalty for each unaligned column. After finding
	a new best overlap, the p-value is calculated by comparing the
	score to the background distribution.

	Targets 0..iq are skipped, and so are their reverse complements, which
	start at len(T_lens) // 2 only when `reverse_complement` is 1.
	"""

	n = len(T_lens) // 2 if reverse_complement == 1 else len(T_lens)
	total_offset = uint64(0)

	for i, nt in enumerate(T_lens):
		nt = uint64(nt)
		results[i, 0] = 1
		results[i, 1] = 0

		if i <= iq or (i >= n and i <= (n + iq)):
			total_offset += nt
			continue

		for k in range(nt+nq-1):
			k = uint64(k)
			t_sums[k] = nq * offset

		for k in range(nt):
			k = uint64(k)
			k_idx = uint64(rr_inv[total_offset + k])
			for l in range(nq):	
				l = uint64(l)
				t_sums[k+l] += gamma[k_idx, l]

		# Only a position holding the maximum can be the final winner: the
		# first one overwrites every field set by an earlier, lower score.
		# Skipping the rest keeps the result and the branches predictable.
		M = t_sums[0]
		for k in range(1, nt+nq-1):
			M = max(M, t_sums[k])

		for k in range(nt+nq-1):
			score = t_sums[k]
			if score != M:
				continue

			overlap = min(k+1, nq) - max(0, k-nt+1)
			if score >= results[i, 1]:
				if score == results[i, 1] and results[i, 2] >= overlap:
					continue

				results[i, 0] = B_cdfs[nt, uint64(score-1)] if score > 0 else 1.0
				results[i, 1] = score
				results[i, 2] = k - nq + 1
				results[i, 3] = overlap

		total_offset += nt


@njit(cache=True)
def _merge_rc_results(results):
	"""An internal method for taking the best across two strands."""

	nt = results.shape[0]
	n = nt // 2
	
	for i in range(n):
		p = min(results[i, 0], results[i+n, 0])
		p = 1 - (1 - p) ** 2

		results[i, 0] = p
		results[i, 4] = 0
		
		if results[i, 1] <= results[i+n, 1]:                
			results[i, 1] = results[i+n, 1]
			results[i, 2] = results[i+n, 2]
			results[i, 3] = results[i+n, 3]
			results[i, 4] = 1
			

@njit(cache=True)
def _merge_rc_results_into(results, out):
	"""Like `_merge_rc_results`, but writes the merged rows into `out`.

	`results` holds both strands; `out` has one row per forward target. Each
	element of `out` is written once, with the same values and tie rule.
	"""

	n = out.shape[0]

	for i in range(n):
		p = min(results[i, 0], results[i+n, 0])
		p = 1 - (1 - p) ** 2

		out[i, 0] = p
		
		if results[i, 1] <= results[i+n, 1]:                
			out[i, 1] = results[i+n, 1]
			out[i, 2] = results[i+n, 2]
			out[i, 3] = results[i+n, 3]
			out[i, 4] = 1
		else:
			out[i, 1] = results[i, 1]
			out[i, 2] = results[i, 2]
			out[i, 3] = results[i, 3]
			out[i, 4] = 0


@njit(parallel=True, cache=True)
def _tomtom(Q, T, Q_lens, T_lens, Q_norm, T_norm, rr_inv, rr_counts, n_nearest, 
	n_score_bins, n_median_bins, n_cache, n_threads, reverse_complement):
	"""An internal function implementing the TOMTOM algorithm.

	This internal function is necessary to handle the numba component of the
	implementation. Here, scratchboard memory is allocated for each thread and
	the main parallel loop is called. Additionally, if reverse complements are
	being considered, values are merged across both strands.
	"""

	T_max = max(T_lens)

	# `_p_values` reads B only at the rows that are target lengths.
	needed = numpy.zeros(T_max+1, dtype=numpy.bool_)
	for t in T_lens:
		needed[t] = True
	
	Q_offsets = numpy.zeros(len(Q_lens)+1, dtype='int64')
	Q_offsets[1:] = numpy.cumsum(Q_lens)
	Q_max = max(Q_lens)
	
	n_in_targets = len(T_lens) // 2 if reverse_complement else len(T_lens)
	n_out_targets = n_in_targets if n_nearest == -1 else n_nearest
	n_outputs = 5 if n_nearest == -1 else 6
	nt = T.shape[-1]

	# Re-usable workspace for each thread instead of re-allocating
	# and freeing large arrays for each example.
	n_len = Q_max*n_score_bins + Q_max*n_cache
	
	_gamma = numpy.empty((n_threads, nt, Q_max), dtype='float64')
	_gamma_int = numpy.empty((n_threads, nt, Q_max), dtype='int16')
	_f = numpy.empty((n_threads, Q_max, n_score_bins+1), dtype='float64')

	# A and A_csum are flat per thread; each query takes a contiguous
	# (nq, nq, n) view so its working set is not strided by n_len.
	_A = numpy.empty((n_threads, Q_max*Q_max*n_len), dtype='float64')
	_B = numpy.empty((n_threads, T_max+1, n_len), dtype='float64')
	_A_csum = numpy.empty((n_threads, Q_max*Q_max*n_len), dtype='float64')

	_medians = numpy.empty((n_threads, Q_max), dtype='float64')
	_median_bins = numpy.empty((n_threads, n_median_bins, 2), dtype='float64')

	_results = numpy.empty((n_threads, len(T_lens), 5), dtype='float64')
	results = numpy.empty((len(Q_lens), n_out_targets, n_outputs), 
		dtype='float64') 

	for i in prange(len(Q_lens)):
		nq = Q_lens[i]
		pid = numba.get_thread_id()

		offset = _integer_distances_and_histogram(Q, T, _gamma[pid], 
			_gamma_int[pid], _f[pid], _medians[pid], _median_bins[pid], Q_norm, 
			T_norm, rr_counts, Q_offsets[i], nq, n_score_bins)

		# The backgrounds span nq*(n_score_bins+offset) bins. When the offset
		# exceeds `n_cache` this can overrun the shared workspace, so allocate
		# a large enough one for this query instead.
		n_needed = nq*n_score_bins + nq*offset
		if n_needed > n_len:
			A = numpy.empty((nq, nq, n_needed), dtype='float64')
			B = numpy.empty((T_max+1, n_needed), dtype='float64')
			A_csum = numpy.empty((nq, nq, n_needed), dtype='float64')
		else:
			n_a = nq*nq*n_needed
			A = _A[pid, :n_a].reshape((nq, nq, n_needed))
			A_csum = _A_csum[pid, :n_a].reshape((nq, nq, n_needed))
			B = _B[pid]

		_p_value_backgrounds(_f[pid], A, B, A_csum, nq, n_score_bins, T_max, 
			offset, needed)

		_p_values(_gamma_int[pid], B, rr_inv, T_lens, -1, nq, offset, 
			_results[pid], reverse_complement)

		# The full-matrix, two-strand case merges straight into the output.
		if reverse_complement == 1 and n_nearest == -1:
			_merge_rc_results_into(_results[pid], results[i])
		else:
			if reverse_complement == 1:
				_merge_rc_results(_results[pid])
			else:
				_results[pid, :, 4] = 0

			if n_nearest == -1:
				results[i] = _results[pid, :n_in_targets]
			else:
				idxs = numpy.argsort(_results[pid, :n_in_targets, 0])[:n_nearest]
				results[i, :, :5] = _results[pid, idxs]
				results[i, :, 5] = idxs


	return results            
  

def tomtom(Qs, Ts, n_nearest=None, n_score_bins=100, n_median_bins=1000, 
	n_target_bins=100, n_cache=100, reverse_complement=True, n_jobs=-1):
	"""A method for assigning p-values to motif similarity.

	This method implements the TOMTOM algorithm for assigning p-values to motif
	similarity scores. TOMTOM accounts for several issues that arise when
	motifs are scanned against each other, including correctly calculating
	scores for overlaps and accounting for motif length and information content 
	within the motifs. 

	At a high level, TOMTOM works by calculating a background distribution of 
	scores for each position in the query and then uses dynamic programming to 
	calculating a distribution of scores for each span of matches, allowing for 
	potential overhangs on either side.

	Importantly, this method implements the "complete score" version of TOMTOM
	which is more robust to edge effects. The "incomplete score" is not a good
	score and so is not implemented. 


	Parameters
	----------
	Qs: list of numpy.ndarray or torch.Tensor with shape (len(alphabet), len)
		A list of query motifs to consider. Each query must have a shape
		according to the PyTorch format where the length is the last aspect.

	Ts: list of numpy.ndarrays or torch.Tensor with shape (len(alphabet), len)
		A list of target motifs to compare each query against. Each target must 
		have a shape according to the PyTorch format where the length is the 
		last aspect.

	n_nearest: int or None, optional
		The number of nearest targets to keep for each query, where nearness is
		defined by the p-value. Setting this can significant reduce memory
		because, otherwise, you get a len(Qs) by len(Ts) complete matrix. If
		None, return the complete matrix. Values larger than len(Ts) are
		clipped to len(Ts). Default is None.

	n_score_bins: int, optional
		The number of bins to use when discretizing scores. A higher number is 
		not necessarily better because you need the data to support each bin
		in the distribution. This is `t` from the TOMTOM paper. Default is 100.

	n_median_bins: int, optional
		The number of bins to use when approximating the medians. More bins
		means higher precision when estimating the median but can also cause it
		to take linearly longer. Default is 1000.

	n_target_bins: int or None, optional
		Whether to use approximate hashing to speed up calculations by merging
		target columns that are similar. This can significantly speed up
		calculations and reduce memory at the cost of approximation. Each value 
		in the columns are binned and targets are merged together if all values 
		fall within the same bins, e.g., if both columns after binning are 
		[5, 11, 0, 1]. This parameter sets the number of bins to use when 
		discretizing the values in the target columns. Fewer bins means more 
		targets get merged together, which can speed up the calculations, but 
		also mean that the resulting p-values are less accurate. Conversely, 
		more bins means that fewer targets get merged together and higher 
		accuracy p-values but slower. If None, don't use approximate hashing.
		Default is 100.

	n_cache: int, optional
		A cache size to use when allocating the scratchpad. A higher number will
		linearly increase the amount of memory used but will not increase the
		amount of compute needed. A query that needs more than this is given
		its own larger scratchpad, so this does not change the results.
		Default is 100.

	reverse_complement: bool, optional
		Whether to automatically compare each query to targets and also the
		reverse complement of the target and merge the scores and p-values
		accordingly. Default is True.

	n_jobs: int, optional
		The number of threads for numba to use when parallelizing the
		processing of query sequences. If -1, use all available threads.
		Default is -1.


	Returns
	-------
	best_p_values: numpy.ndarray, shape=(len(Qs), len(Ts))
		The p-value of the best alignment between each query and each target.

	best_scores: numpy.ndarray, shape=(len(Qs), len(Ts))
		The scores of the best alignment between each query and each target.

	best_offsets: numpy.ndarray, shape=(len(Qs), len(Ts))
		The offset of the best alignment between each query and each target.

	best_overlaps: numpy.ndarray, shape=(len(Qs), len(Ts))
		The overlap of the best alignment between each query and each target.

	best_strands: numpy.ndarray, shape=(len(Qs), len(Ts))
		The strand for the best alignment between each query and each target.

	best_idxs: numpy.ndarray, shape=(len(Qs), len(Ts)), optional
		When returning only a number of nearest neighbors, the index in the
		original ordering of the targets corresponding to each returned
		neighbor. These will be sorted by p-value.
	"""

	if n_jobs != -1:
		_n_jobs = numba.get_num_threads()
		numba.set_num_threads(n_jobs)
	else:
		n_jobs = _n_jobs = numba.get_num_threads()

	# Asking for more neighbors than there are targets returns every target;
	# the surplus columns would otherwise be left uninitialized.
	if n_nearest is None:
		n_nearest = -1
	else:
		n_nearest = min(n_nearest, len(Ts))

	if not isinstance(Qs[0], numpy.ndarray):
		Qs = [Q.numpy() for Q in Qs]
	
	if not isinstance(Ts[0], numpy.ndarray):
		Ts = [T.numpy() for T in Ts]

	Q_lens = numpy.array([Q.shape[-1] for Q in Qs], dtype='int64')
	Q = numpy.concatenate(Qs, axis=-1)
	Q_norm = (Q ** 2).sum(axis=0)
	
	if reverse_complement:        
		Ts = Ts + [T[::-1, ::-1] for T in Ts]
	
	T_lens = numpy.array([T.shape[-1] for T in Ts], dtype='int64')
	T = numpy.concatenate(Ts, axis=-1)
	T_norm = (T ** 2).sum(axis=0)

	if Q_norm.max() == 0 or T_norm.max() == 0:
		raise ValueError("Cannot have all-zeroes as targets or query.")

	if n_target_bins is not None:
		T_min = T.min(axis=-1, keepdims=True)
		T_max = T.max(axis=-1, keepdims=True)
		T_max[T_max == T_min] = T_min[T_max == T_min] + 1

		T_ints = numpy.around((T - T_min) / (T_max - T_min) * (n_target_bins-1))
		T_ints = T_ints.T.dot(n_target_bins ** numpy.arange(len(T))[:, None])
		_, rr_idxs, rr_inv, rr_counts = numpy.unique(T_ints.flatten(), 
			return_index=True, return_inverse=True, return_counts=True)

		T = T[:, rr_idxs]
		T_norm = T_norm[rr_idxs]
		rr_inv = rr_inv.astype('uint64')
	else:
		rr_inv = numpy.arange(T.shape[-1])
		rr_counts = numpy.ones_like(rr_inv)
	
	###
	
	results = _tomtom(Q, T, Q_lens, T_lens, Q_norm, T_norm, rr_inv, rr_counts, 
		n_nearest, n_score_bins, n_median_bins, n_cache, n_jobs, 
		int(reverse_complement))

	if n_jobs != -1:
		numba.set_num_threads(_n_jobs)

	return results.transpose(2, 0, 1)
