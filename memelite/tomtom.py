# tomtom.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com> 

import time
import math
import numpy
import numba

from numba import njit
from numba import prange
from numba.core.cpu_options import ParallelOptions
from numpy import uint64
from numpy import int64
from numpy import int32


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

	# Every value is x_min when the range is zero, so x_min is the median.
	if x_max == 0:
		return x_min

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


# Bit position of a one-bit uint64 v: _DEBRUIJN[(v * 0x03f79d71b4cb0a89) >> 58].
_DEBRUIJN = numpy.array([0, 1, 48, 2, 57, 49, 28, 3, 61, 58, 50, 42, 38, 29,
	17, 4, 62, 55, 59, 36, 53, 51, 43, 22, 45, 39, 33, 30, 24, 18, 12, 5, 63, 47,
	56, 27, 60, 41, 37, 16, 54, 35, 52, 21, 44, 32, 23, 11, 46, 26, 40, 15, 34,
	20, 31, 10, 25, 14, 19, 9, 13, 8, 7, 6], dtype=numpy.int64)


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
	a check that vectorizes (a min over zb ^ b, zero only on a hit). In a hit
	block the elements in the bin are read from a 64-bit mask, lowest bit
	first, so they are visited in ascending order without a data-dependent
	branch per element. Bitwise equal to `_binned_median`. `bins` must be
	C-contiguous.
	"""

	n, n_bins = len(x), len(bins)
	x_max -= x_min

	# Every value is x_min when the range is zero, so x_min is the median.
	if x_max == 0:
		return x_min

	for i in range(n):
		zb[i] = int((x[i] - x_min) / x_max * (n_bins - 1))

	cnt = bins.reshape(-1).view(numpy.int64)[:n_bins]
	cnt[:] = 0
	# Every index lies in [0, n_bins), so an unsigned index is the same bin
	# and drops numba's negative-index wraparound from the scatter.
	for i in range(n):
		cnt[uint64(zb[i])] += counts[i]

	count = 0
	for b in range(n_bins):
		count += cnt[b]
		if count >= halfway:
			s = 0.0
			n_64 = n - n % 64
			b32 = numpy.int32(b)
			for i0 in range(0, n_64, 64):
				d = numpy.uint32(0xFFFFFFFF)
				for i in range(i0, i0 + 64):
					d = min(d, numpy.uint32(zb[i] ^ b32))
				if d != 0:
					continue

				m = uint64(0)
				for t in range(64):
					m |= uint64(zb[i0 + t] == b32) << uint64(t)
				while m != 0:
					low = m & (~m + uint64(1))
					i = i0 + _DEBRUIJN[(low * uint64(0x03f79d71b4cb0a89)) >> uint64(58)]
					s += x[i] * counts[i]
					m ^= low

			for i in range(n_64, n):
				if zb[i] == b:
					s += x[i] * counts[i]

			return s / cnt[b]

	return -99999


@njit(cache=True)
def _median_from_counts(x, cnt, zb, counts, halfway):
	"""The scan and stopping-bin sum of `_binned_median_z`, given its counts
	`cnt` and bin indices `zb`: the same operations in the same order."""

	n, n_bins = len(x), len(cnt)
	count = 0
	for b in range(n_bins):
		count += cnt[b]
		if count >= halfway:
			s = 0.0
			n_64 = n - n % 64
			b32 = numpy.int32(b)
			for i0 in range(0, n_64, 64):
				d = numpy.uint32(0xFFFFFFFF)
				for i in range(i0, i0 + 64):
					d = min(d, numpy.uint32(zb[i] ^ b32))
				if d != 0:
					continue

				m = uint64(0)
				for t in range(64):
					m |= uint64(zb[i0 + t] == b32) << uint64(t)
				while m != 0:
					low = m & (~m + uint64(1))
					i = i0 + _DEBRUIJN[(low * uint64(0x03f79d71b4cb0a89)) >> uint64(58)]
					s += x[i] * counts[i]
					m ^= low

			for i in range(n_64, n):
				if zb[i] == b:
					s += x[i] * counts[i]

			return s / cnt[b]

	return -99999


@njit(cache=True)
def _binned_median_block4(r0, r1, r2, r3, mn0, mx0, mn1, mx1, mn2, mx2, mn3,
	mx3, bins, counts, counts32, zb4, halfway):
	"""`_binned_median_z` of four rows, bitwise the same.

	Each row's bin indices are `_binned_median_z`'s expression, in its own
	loop. The four count histograms are filled in one loop over the targets,
	so each target's count is loaded once for four scatters. They count in
	int32 from `counts32`, the same counts; the caller passes `counts32`
	only when sum(counts) fits int32, so no bin or partial count overflows.
	The histograms live in `bins`'s own memory, which holds 4 x n_bins
	int32. The scan and the stopping-bin sum are `_binned_median_z`'s, per
	row, and read the int64 `counts`.
	"""

	n, n_bins = len(r0), len(bins)
	z0, z1, z2, z3 = zb4[0], zb4[1], zb4[2], zb4[3]

	# A row whose values are all equal has a zero range. Its elements all go
	# in bin 0, and its median is its minimum, as in `_binned_median_z`.
	mx0 -= mn0
	if mx0 == 0:
		z0[:n] = 0
	else:
		for i in range(n):
			z0[i] = int((r0[i] - mn0) / mx0 * (n_bins - 1))
	mx1 -= mn1
	if mx1 == 0:
		z1[:n] = 0
	else:
		for i in range(n):
			z1[i] = int((r1[i] - mn1) / mx1 * (n_bins - 1))
	mx2 -= mn2
	if mx2 == 0:
		z2[:n] = 0
	else:
		for i in range(n):
			z2[i] = int((r2[i] - mn2) / mx2 * (n_bins - 1))
	mx3 -= mn3
	if mx3 == 0:
		z3[:n] = 0
	else:
		for i in range(n):
			z3[i] = int((r3[i] - mn3) / mx3 * (n_bins - 1))

	c4 = bins.reshape(-1).view(numpy.int32)
	c0 = c4[0:n_bins]
	c1 = c4[n_bins:2 * n_bins]
	c2 = c4[2 * n_bins:3 * n_bins]
	c3 = c4[3 * n_bins:4 * n_bins]
	for b in range(4 * n_bins):
		c4[b] = 0
	# Every index lies in [0, n_bins); see `_binned_median_z`.
	for i in range(n):
		ci = counts32[i]
		c0[uint64(z0[i])] += ci
		c1[uint64(z1[i])] += ci
		c2[uint64(z2[i])] += ci
		c3[uint64(z3[i])] += ci

	m0 = _median_from_counts(r0, c0, z0, counts, halfway)
	m1 = _median_from_counts(r1, c1, z1, counts, halfway)
	m2 = _median_from_counts(r2, c2, z2, counts, halfway)
	m3 = _median_from_counts(r3, c3, z3, counts, halfway)
	return (mn0 if mx0 == 0 else m0, mn1 if mx1 == 0 else m1,
		mn2 if mx2 == 0 else m2, mn3 if mx3 == 0 else m3)


@njit(cache=True, inline='always')
def _column_distances(X, c, Y, Y_norm, xn, g_row, x2, mxb, mnb, zb, 
	median_bins, Y_counts, halfway):
	"""One query column's distance row, its min and max, and its median.

	Writes -sqrt(distance) between query column `c` of `X` and every column
	of `Y` into `g_row` and returns (min, max, binned median) of that row.
	Everything depends only on the column's values and the targets, so a
	column that occurs more than once gives the same row and scalars.

	The distance terms are subtracted in the same order as before; `2 * X` is
	hoisted, which is the same product. The min/max pass is separate so the
	distance loop vectorizes. It runs as 64 lanes through numpy.maximum/
	minimum, which compile to packed code where a serial max/min chain
	cannot. Every g value is -sqrt(z) with z > 0 or +0.0: never NaN and never
	-0.0, so max and min give the same bits in any order.
	"""

	n_a, n_y = Y.shape[0], Y.shape[-1]

	if n_a == 4:
		x0 = 2 * X[0, c]
		x1 = 2 * X[1, c]
		x2_ = 2 * X[2, c]
		x3 = 2 * X[3, c]
		for j in range(n_y):
			z = xn + Y_norm[j]
			z -= x0 * Y[0, j]
			z -= x1 * Y[1, j]
			z -= x2_ * Y[2, j]
			z -= x3 * Y[3, j]
			g_row[j] = -math.sqrt(z) if z > 0 else 0
	else:
		for k in range(n_a):
			x2[k] = 2 * X[k, c]
		for j in range(n_y):
			z = xn + Y_norm[j]
			for k in range(n_a):
				z -= x2[k] * Y[k, j]
			g_row[j] = -math.sqrt(z) if z > 0 else 0

	return _column_stats(g_row, mxb, mnb, zb, median_bins, Y_counts, halfway)


@njit(cache=True, inline='always')
def _column_stats(g_row, mxb, mnb, zb, median_bins, Y_counts, halfway):
	"""(min, max, binned median) of one distance row; see `_column_distances`."""

	n_y = len(g_row)
	z_min_, z_max_ = 9999999.9, -9999999.9
	nb = n_y - n_y % 64
	if nb > 0:
		for u in range(64):
			mxb[u] = g_row[u]
			mnb[u] = g_row[u]
		for j in range(64, nb, 64):
			numpy.maximum(mxb, g_row[j:j+64], mxb)
			numpy.minimum(mnb, g_row[j:j+64], mnb)
		for u in range(64):
			z_max_ = max(z_max_, mxb[u])
			z_min_ = min(z_min_, mnb[u])
	for j in range(nb, n_y):
		z = g_row[j]
		z_max_ = max(z_max_, z)
		z_min_ = min(z_min_, z)

	m = _binned_median_z(g_row, median_bins, z_min_, z_max_, Y_counts,
		zb, halfway)
	return z_min_, z_max_, m


@njit(cache=True, inline='always')
def _distance4(xn, x0, x1, x2, x3, yn, y0, y1, y2, y3):
	"""One element of a distance row: the 4-letter loop's operations, in its
	order, from `_column_distances`."""

	z = xn + yn
	z -= x0 * y0
	z -= x1 * y1
	z -= x2 * y2
	z -= x3 * y3
	return -math.sqrt(z) if z > 0 else 0.0


@njit(cache=True)
def _distances_block4(X, Y, Y_norm, X_norm, c0, c1, c2, c3, r0, r1, r2, r3,
	L):
	"""Distance rows of four query columns in one sweep over the targets,
	with each row's min and max.

	Each target column's four values and norm are loaded once for the four
	query columns. Every row element gets exactly the operations of the
	4-letter loop in `_column_distances`, in the same order, so the rows are
	bitwise the same; only the loop nest is reorganized.

	Returns (min0, max0, min1, max1, min2, max2, min3, max3). The targets
	run in tiles of 64, and element t of a tile updates lane t of each row's
	running max and min in `L` (8 x 64 scratch) in the same loop that writes
	it, which vectorizes: the lanes are 8 masked read-modify-writes per
	8 targets that stay in L1. A separate pass over each row (as
	`_column_stats` does) re-reads 4 rows of 150 kB. Every value is -sqrt(z)
	with z > 0 or +0.0, never NaN and never -0.0, so max and min over the
	lanes, the tail and the +-9999999.9 starting values give the bits
	`_column_stats` gives, in any order.

	Y must have four rows. For an F-ordered Y (`symmetric_tomtom`, or
	`tomtom` with n_target_bins=None) the check below is what fixes the
	column stride at four values; without it the stride is a runtime value,
	the loop stays scalar, and a whole tomtom call was 0.37 s slower.
	"""

	n_a, n_y = Y.shape[0], Y.shape[-1]
	lo, hi = 9999999.9, -9999999.9
	if n_a != 4:
		return lo, hi, lo, hi, lo, hi, lo, hi

	a0, a1, a2, a3 = 2 * X[0, c0], 2 * X[1, c0], 2 * X[2, c0], 2 * X[3, c0]
	b0, b1, b2, b3 = 2 * X[0, c1], 2 * X[1, c1], 2 * X[2, c1], 2 * X[3, c1]
	d0, d1, d2, d3 = 2 * X[0, c2], 2 * X[1, c2], 2 * X[2, c2], 2 * X[3, c2]
	e0, e1, e2, e3 = 2 * X[0, c3], 2 * X[1, c3], 2 * X[2, c3], 2 * X[3, c3]
	xa, xb, xd, xe = X_norm[c0], X_norm[c1], X_norm[c2], X_norm[c3]

	M0, N0, M1, N1, M2, N2, M3, N3 = L[0], L[1], L[2], L[3], L[4], L[5], \
		L[6], L[7]
	for t in range(64):
		M0[t], M1[t], M2[t], M3[t] = hi, hi, hi, hi
		N0[t], N1[t], N2[t], N3[t] = lo, lo, lo, lo

	nb = n_y - n_y % 64
	for j0 in range(0, nb, 64):
		for t in range(64):
			j = j0 + t
			y0, y1, y2, y3, yn = Y[0, j], Y[1, j], Y[2, j], Y[3, j], Y_norm[j]

			g = _distance4(xa, a0, a1, a2, a3, yn, y0, y1, y2, y3)
			r0[j] = g
			M0[t] = g if g > M0[t] else M0[t]
			N0[t] = g if g < N0[t] else N0[t]

			g = _distance4(xb, b0, b1, b2, b3, yn, y0, y1, y2, y3)
			r1[j] = g
			M1[t] = g if g > M1[t] else M1[t]
			N1[t] = g if g < N1[t] else N1[t]

			g = _distance4(xd, d0, d1, d2, d3, yn, y0, y1, y2, y3)
			r2[j] = g
			M2[t] = g if g > M2[t] else M2[t]
			N2[t] = g if g < N2[t] else N2[t]

			g = _distance4(xe, e0, e1, e2, e3, yn, y0, y1, y2, y3)
			r3[j] = g
			M3[t] = g if g > M3[t] else M3[t]
			N3[t] = g if g < N3[t] else N3[t]

	mn0, mx0, mn1, mx1, mn2, mx2, mn3, mx3 = lo, hi, lo, hi, lo, hi, lo, hi
	for j in range(nb, n_y):
		y0, y1, y2, y3, yn = Y[0, j], Y[1, j], Y[2, j], Y[3, j], Y_norm[j]
		g = _distance4(xa, a0, a1, a2, a3, yn, y0, y1, y2, y3)
		r0[j] = g
		mn0, mx0 = min(mn0, g), max(mx0, g)
		g = _distance4(xb, b0, b1, b2, b3, yn, y0, y1, y2, y3)
		r1[j] = g
		mn1, mx1 = min(mn1, g), max(mx1, g)
		g = _distance4(xd, d0, d1, d2, d3, yn, y0, y1, y2, y3)
		r2[j] = g
		mn2, mx2 = min(mn2, g), max(mx2, g)
		g = _distance4(xe, e0, e1, e2, e3, yn, y0, y1, y2, y3)
		r3[j] = g
		mn3, mx3 = min(mn3, g), max(mx3, g)

	for t in range(64):
		mx0, mn0 = max(mx0, M0[t]), min(mn0, N0[t])
		mx1, mn1 = max(mx1, M1[t]), min(mn1, N1[t])
		mx2, mn2 = max(mx2, M2[t]), min(mn2, N2[t])
		mx3, mn3 = max(mx3, M3[t]), min(mn3, N3[t])

	return mn0, mx0, mn1, mx1, mn2, mx2, mn3, mx3


@njit(cache=True)
def _distances_block2(X, Y, Y_norm, X_norm, c0, c1, r0, r1, L):
	"""`_distances_block4` for two query columns, used for the columns left
	over: returns (min0, max0, min1, max1), with each row's running max and
	min kept in lanes of `L` in the loop that writes it, as there. A single
	column passes c1 = c0 and a spare row as r1."""

	n_a, n_y = Y.shape[0], Y.shape[-1]
	lo, hi = 9999999.9, -9999999.9
	if n_a != 4:
		return lo, hi, lo, hi

	a0, a1, a2, a3 = 2 * X[0, c0], 2 * X[1, c0], 2 * X[2, c0], 2 * X[3, c0]
	b0, b1, b2, b3 = 2 * X[0, c1], 2 * X[1, c1], 2 * X[2, c1], 2 * X[3, c1]
	xa, xb = X_norm[c0], X_norm[c1]

	M0, N0, M1, N1 = L[0], L[1], L[2], L[3]
	for t in range(64):
		M0[t], M1[t] = hi, hi
		N0[t], N1[t] = lo, lo

	nb = n_y - n_y % 64
	for j0 in range(0, nb, 64):
		for t in range(64):
			j = j0 + t
			y0, y1, y2, y3, yn = Y[0, j], Y[1, j], Y[2, j], Y[3, j], Y_norm[j]

			g = _distance4(xa, a0, a1, a2, a3, yn, y0, y1, y2, y3)
			r0[j] = g
			M0[t] = g if g > M0[t] else M0[t]
			N0[t] = g if g < N0[t] else N0[t]

			g = _distance4(xb, b0, b1, b2, b3, yn, y0, y1, y2, y3)
			r1[j] = g
			M1[t] = g if g > M1[t] else M1[t]
			N1[t] = g if g < N1[t] else N1[t]

	mn0, mx0, mn1, mx1 = lo, hi, lo, hi
	for j in range(nb, n_y):
		y0, y1, y2, y3, yn = Y[0, j], Y[1, j], Y[2, j], Y[3, j], Y_norm[j]
		g = _distance4(xa, a0, a1, a2, a3, yn, y0, y1, y2, y3)
		r0[j] = g
		mn0, mx0 = min(mn0, g), max(mx0, g)
		g = _distance4(xb, b0, b1, b2, b3, yn, y0, y1, y2, y3)
		r1[j] = g
		mn1, mx1 = min(mn1, g), max(mx1, g)

	for t in range(64):
		mx0, mn0 = max(mx0, M0[t]), min(mn0, N0[t])
		mx1, mn1 = max(mx1, M1[t]), min(mn1, N1[t])

	return mn0, mx0, mn1, mx1


@njit(cache=True)
def _halfway(Y_counts):
	halfway = 0
	for j in range(len(Y_counts)):
		halfway += Y_counts[j]
	return halfway / 2


@njit(cache=True)
def _fill_column_cache(X, Y, X_norm, Y_norm, Y_counts, cols, G_cache, 
	S_cache, n_median_bins):
	"""Distance rows and (min, max, median) for the columns in `cols`.

	Row s of `G_cache` and `S_cache` hold column `cols[s]`'s results from
	`_column_distances`, the same function the per-query path runs, so a
	query can read them in place of recomputing its own. Filled once, before
	the parallel loop, and only read afterwards.
	"""

	n_a, n_y = Y.shape[0], Y.shape[-1]
	x2 = numpy.empty(n_a, dtype=numpy.float64)
	zb = numpy.empty(n_y, dtype=numpy.int32)
	mxb = numpy.empty(64, dtype=numpy.float64)
	mnb = numpy.empty(64, dtype=numpy.float64)
	median_bins = numpy.empty((n_median_bins, 2), dtype=numpy.float64)
	halfway = _halfway(Y_counts)

	for s in range(len(cols)):
		c = cols[s]
		z_min_, z_max_, m = _column_distances(X, c, Y, Y_norm, X_norm[c], 
			G_cache[s], x2, mxb, mnb, zb, median_bins, Y_counts, halfway)
		S_cache[s, 0] = z_min_
		S_cache[s, 1] = z_max_
		S_cache[s, 2] = m


@njit(cache=True)
def _binned_load(h_int, h_f, gamma_int, f_row, k):
	"""Copy a cached `gamma_int` column and `f` row into this query's."""

	for j in range(len(h_int)):
		gamma_int[j, k] = h_int[j]
	for x in range(len(h_f)):
		f_row[x] = h_f[x]


@njit(cache=True)
def _binned_save(h_int, h_f, gamma_int, f_row, k):
	"""Save a `gamma_int` column and `f` row into the cache."""

	for j in range(len(h_int)):
		h_int[j] = gamma_int[j, k]
	for x in range(len(h_f)):
		h_f[x] = f_row[x]


@njit(cache=True, inline='always')
def _binned_column(row, mi, bin_scale, offset, w, zb, gamma_int, k, f_row):
	"""The binned stage of one query column: its bin indices, its
	`gamma_int[:, k]` column and its histogram row `f_row`."""

	n_y = len(w)
	for j in range(n_y):
		zb[j] = math.floor((row[j] - mi) * bin_scale + 0.5)

	for j in range(n_y):
		x = zb[j]
		f_row[uint64(x)] += w[j]
		gamma_int[j, k] = x - offset


@njit(cache=True)
def _binned_block4(r0, r1, r2, r3, m0, m1, m2, m3, bin_scale, offset, w,
	zb4, gamma_int, k0, k1, k2, k3, f0, f1, f2, f3):
	"""The binned stage of four query columns in one sweep over the targets.

	Each row's bin indices are the expression of the single-column loop in
	`_integer_distances_and_histogram`. The store/scatter loop then handles
	the four columns per target: their `gamma_int[j, k]` stores share one
	row, and each `f` row receives its additions in ascending j, as before,
	so every sum is bitwise the same.

	Each column is finished (scatter, then store) before the next is
	loaded. With the four loads, four stores and four scatters grouped
	instead, LLVM kept all four bin indices live and reloaded pointers from
	the stack 9 times per target; this order holds one index at a time.
	`_binned_block4c` is the faster form for consecutive columns.
	"""

	n_y = len(w)
	z0, z1, z2, z3 = zb4[0], zb4[1], zb4[2], zb4[3]
	for j in range(n_y):
		z0[j] = math.floor((r0[j] - m0) * bin_scale + 0.5)
	for j in range(n_y):
		z1[j] = math.floor((r1[j] - m1) * bin_scale + 0.5)
	for j in range(n_y):
		z2[j] = math.floor((r2[j] - m2) * bin_scale + 0.5)
	for j in range(n_y):
		z3[j] = math.floor((r3[j] - m3) * bin_scale + 0.5)

	for j in range(n_y):
		wj = w[j]
		x = z0[j]
		f0[uint64(x)] += wj
		gamma_int[j, k0] = x - offset
		x = z1[j]
		f1[uint64(x)] += wj
		gamma_int[j, k1] = x - offset
		x = z2[j]
		f2[uint64(x)] += wj
		gamma_int[j, k2] = x - offset
		x = z3[j]
		f3[uint64(x)] += wj
		gamma_int[j, k3] = x - offset


@njit(cache=True)
def _binned_block4c(r0, r1, r2, r3, m0, m1, m2, m3, bin_scale, offset, w,
	zb4, gamma_int, kb, f0, f1, f2, f3):
	"""`_binned_block4` for four consecutive query columns, whose `gamma_int`
	columns are kb + 3, kb + 2, kb + 1 and kb. One unsigned base column lets
	the four stores share one row pointer at constant displacements, where
	four independent k need four offsets and spill. The caller must decide
	consecutiveness: an equality test on k0..k3 inside one kernel lets LLVM
	rewrite kb + 3 back into k0 and the four offsets return."""

	n_y = len(w)
	z0, z1, z2, z3 = zb4[0], zb4[1], zb4[2], zb4[3]
	for j in range(n_y):
		z0[j] = math.floor((r0[j] - m0) * bin_scale + 0.5)
	for j in range(n_y):
		z1[j] = math.floor((r1[j] - m1) * bin_scale + 0.5)
	for j in range(n_y):
		z2[j] = math.floor((r2[j] - m2) * bin_scale + 0.5)
	for j in range(n_y):
		z3[j] = math.floor((r3[j] - m3) * bin_scale + 0.5)

	k0, k1, k2, k3 = kb + uint64(3), kb + uint64(2), kb + uint64(1), kb

	for j in range(n_y):
		wj = w[j]
		x = z0[j]
		f0[uint64(x)] += wj
		gamma_int[j, k0] = x - offset
		x = z1[j]
		f1[uint64(x)] += wj
		gamma_int[j, k1] = x - offset
		x = z2[j]
		f2[uint64(x)] += wj
		gamma_int[j, k2] = x - offset
		x = z3[j]
		f3[uint64(x)] += wj
		gamma_int[j, k3] = x - offset


@njit(cache=True)
def _binned_block2(r0, r1, m0, m1, bin_scale, offset, w,
	zb4, gamma_int, k0, k1, f0, f1):
	"""`_binned_block4` for two query columns, used for a leftover pair."""

	n_y = len(w)
	z0, z1 = zb4[0], zb4[1]
	for j in range(n_y):
		z0[j] = math.floor((r0[j] - m0) * bin_scale + 0.5)
	for j in range(n_y):
		z1[j] = math.floor((r1[j] - m1) * bin_scale + 0.5)
	for j in range(n_y):
		wj = w[j]
		x = z0[j]
		f0[uint64(x)] += wj
		gamma_int[j, k0] = x - offset
		x = z1[j]
		f1[uint64(x)] += wj
		gamma_int[j, k1] = x - offset


@njit(cache=True, inline='always')
def _fits_gamma(offset, n_bins, g_max):
	"""Whether x - offset lies in [-g_max-1, g_max] for every x in
	[0, n_bins]. `offset` is the signed value (int64 of the returned uint64)."""

	return int64(offset) <= int64(g_max) + 1 and \
		int64(n_bins) - int64(offset) <= int64(g_max)


@njit(cache=True)
def _integer_distances_and_histogram(X, Y, gamma, gamma_int, f, medians, 
	median_bins, X_norm, Y_norm, Y_counts, nq_csum, nq, n_bins, q_slot=None,
	G_cache=None, S_cache=None, H_keys=None, H_filled=None, H_int=None,
	H_f=None, w_pre=None, halfway_pre=None, g_max=-1):
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

	When `q_slot` is given, a query column c with q_slot[c] >= 0 takes its
	distance row, min, max and median from row q_slot[c] of `G_cache` and
	`S_cache` (see `_fill_column_cache`) instead of computing them. Those are
	the values the computation would give, bit for bit.

	When `H_keys` is also given, the binned stage of the cached columns with
	slot < H_int.shape[1] is cached too. A column's `gamma_int` column and
	`f` row depend only on its slot and the query's (i_min, bin_scale), which
	also fix `offset`. H_keys[0, 0] counts the (i_min, bin_scale) classes
	seen, stored in H_keys[1:]. The first query of a class to reach a slot
	computes both as usual and saves them to H_int[h, slot] and H_f[h, slot];
	later ones copy them. The arrays belong to one thread.

	`w_pre` and `halfway_pre`, when given, are the target weights
	Y_counts / sum(Y_counts) and `_halfway(Y_counts)`, computed once by the
	caller with the same expressions used here; they depend only on the
	targets, so every query reads the same values.

	`g_max` >= 0 says `gamma_int` holds only [-g_max-1, g_max] (127 for
	int8). Every stored value is x - offset with x in [0, n_bins] (x also
	indexes `f`, which has n_bins + 1 columns), so they all fit when
	`_fits_gamma(offset, n_bins, g_max)`. When they do not, the function
	returns the offset right after computing it, writing neither `gamma_int`
	nor `f`, and the caller reruns it into a wider `gamma_int`.
	"""

	i_min, bin_scale, offset = _distances_and_medians(X, Y, gamma, medians,
		median_bins, X_norm, Y_norm, Y_counts, nq_csum, nq, n_bins, q_slot,
		S_cache, halfway_pre)
	if g_max >= 0 and not _fits_gamma(offset, n_bins, g_max):
		return uint64(offset)

	return _integer_histogram(Y, gamma, gamma_int, f, medians, Y_counts,
		nq_csum, nq, i_min, bin_scale, offset, q_slot, G_cache, H_keys,
		H_filled, H_int, H_f, w_pre)


@njit(cache=True)
def _distances_and_medians(X, Y, gamma, medians, median_bins, X_norm, Y_norm,
	Y_counts, nq_csum, nq, n_bins, q_slot=None, S_cache=None,
	halfway_pre=None, Y_counts32=None):
	"""Steps (1) and (2) of `_integer_distances_and_histogram`: each query
	column's distance row (into `gamma`) and median (into `medians`), and
	the query's (i_min, bin_scale, offset).

	None of it touches `gamma_int`, so it is compiled once, while
	`_integer_histogram` is compiled for each `gamma_int` type (int8 batched,
	int16 alone). A query that does not fit the int8 batch keeps these
	results and reruns only `_integer_histogram`.
	"""

	# `gamma` is private scratch, so its buffer is used as (n_rows, n_y): each
	# query column's distances are then contiguous for every loop below.
	n_a, n_y = Y.shape[0], Y.shape[-1]
	g = gamma.reshape((gamma.shape[1], gamma.shape[0]))
	x2 = numpy.empty(n_a, dtype=numpy.float64)
	zb4 = numpy.empty((4, n_y), dtype=numpy.int32)
	zb = zb4[0]
	mxb = numpy.empty(64, dtype=numpy.float64)
	mnb = numpy.empty(64, dtype=numpy.float64)
	L = numpy.empty((8, 64), dtype=numpy.float64)
	if halfway_pre is None:
		halfway = _halfway(Y_counts)
	else:
		halfway = halfway_pre

	# Each column's (min, max, median) goes into smin/smax/medians, and the
	# reduction over columns runs afterwards in column order. With four
	# letters, columns without a cached row are computed four at a time by
	# `_distances_block4`, so one sweep over the targets serves four query
	# columns and gives their min and max; each block's rows get their median
	# right away, while they are still in cache. Of the up to three columns
	# left over, a pair goes through `_distances_block2`, which also gives
	# their min and max, and a single one goes through it paired with itself.
	#
	# With an alphabet other than four letters, uncached columns are computed
	# one at a time: they are collected in `rest` and share one call of the
	# inlined `_column_distances`, so its code is compiled once rather than
	# twice. Each column's results depend only on that column, so the order
	# they are computed in does not matter.
	smin = numpy.empty(nq, dtype=numpy.float64)
	smax = numpy.empty(nq, dtype=numpy.float64)
	todo = numpy.empty(4, dtype=numpy.int64)
	rest = numpy.empty(nq, dtype=numpy.int64)
	n_todo = 0
	n_rest = 0
	for i in range(nq):
		c = i + nq_csum
		slot = -1
		if q_slot is not None:
			slot = q_slot[c]

		# `S_cache is not None` first, so numba removes the branch while typing
		# when there is no cache; numba 0.60 does not infer that from `slot`.
		if S_cache is not None and slot >= 0:
			smin[i] = S_cache[slot, 0]
			smax[i] = S_cache[slot, 1]
			medians[i] = S_cache[slot, 2]
		elif n_a == 4:
			todo[n_todo] = i
			n_todo += 1
			if n_todo == 4:
				i0, i1, i2, i3 = todo[0], todo[1], todo[2], todo[3]
				(smin[i0], smax[i0], smin[i1], smax[i1], smin[i2], smax[i2],
					smin[i3], smax[i3]) = _distances_block4(X, Y, Y_norm,
					X_norm, i0 + nq_csum, i1 + nq_csum, i2 + nq_csum,
					i3 + nq_csum, g[i0], g[i1], g[i2], g[i3], L)
				if Y_counts32 is not None and len(Y_counts32) == n_y:
					(medians[i0], medians[i1], medians[i2],
						medians[i3]) = _binned_median_block4(g[i0], g[i1],
						g[i2], g[i3], smin[i0], smax[i0], smin[i1], smax[i1],
						smin[i2], smax[i2], smin[i3], smax[i3], median_bins,
						Y_counts, Y_counts32, zb4, halfway)
				else:
					for u in range(4):
						iu = todo[u]
						medians[iu] = _binned_median_z(g[iu], median_bins, 
							smin[iu], smax[iu], Y_counts, zb, halfway)
				n_todo = 0
		else:
			rest[n_rest] = i
			n_rest += 1

	if n_todo >= 2:
		i0, i1 = todo[0], todo[1]
		smin[i0], smax[i0], smin[i1], smax[i1] = _distances_block2(X, Y,
			Y_norm, X_norm, i0 + nq_csum, i1 + nq_csum, g[i0], g[i1], L)
		for u in range(2):
			iu = todo[u]
			medians[iu] = _binned_median_z(g[iu], median_bins,
				smin[iu], smax[iu], Y_counts, zb, halfway)
		todo[0] = todo[2]
		n_todo -= 2
	if n_todo == 1:
		# A single column runs as a pair with itself, the second row going
		# to a spare row whose results are discarded.
		i0 = todo[0]
		if nq < g.shape[0]:
			r1 = g[nq]
		else:
			r1 = numpy.empty(n_y, dtype=numpy.float64)
		smin[i0], smax[i0], _mn, _mx = _distances_block2(X, Y, Y_norm, X_norm,
			i0 + nq_csum, i0 + nq_csum, g[i0], r1, L)
		medians[i0] = _binned_median_z(g[i0], median_bins, smin[i0],
			smax[i0], Y_counts, zb, halfway)
		n_todo = 0
	for u in range(n_rest):
		iu = rest[u]
		c = iu + nq_csum
		smin[iu], smax[iu], medians[iu] = _column_distances(X, c, Y, Y_norm,
			X_norm[c], g[iu], x2, mxb, mnb, zb, median_bins, Y_counts, halfway)

	z_min, z_max = 9999999.9, -9999999.9
	for i in range(nq):
		m = medians[i]
		z_min = min(z_min, smin[i] - m)
		z_max = max(z_max, smax[i] - m)
			
	# Find the minimum value and the number of bins needed to get there.
	# z_max - i_min is below 1 only when z_min is 0, i.e. every column's
	# median is its minimum, and z_max < 1. When every target column scores
	# the same it is 0 or round-off, and dividing by it would stretch that
	# round-off over all n_bins bins, so the divisor is at least 1.
	i_min = int(math.floor(z_min)) #offset
	bin_scale = int(math.floor(n_bins / max(z_max - i_min, 1.0))) #scale
	offset = -i_min * bin_scale
	return i_min, bin_scale, offset


@njit(cache=True)
def _integer_histogram(Y, gamma, gamma_int, f, medians, Y_counts, nq_csum, nq,
	i_min, bin_scale, offset, q_slot=None, G_cache=None, H_keys=None,
	H_filled=None, H_int=None, H_f=None, w_pre=None):
	"""Steps (3) and (4) of `_integer_distances_and_histogram`, from the
	outputs of `_distances_and_medians`: `gamma_int` and the histogram `f`.
	Adds i_min to `medians`, so it runs once per `_distances_and_medians`.
	"""

	n_a, n_y = Y.shape[0], Y.shape[-1]
	g = gamma.reshape((gamma.shape[1], gamma.shape[0]))
	todo = numpy.empty(4, dtype=numpy.int64)

	for i in range(nq):
		medians[i] = medians[i] + i_min
	
	f[:] = 0
	if w_pre is None:
		ys = numpy.sum(Y_counts)
		w = numpy.empty(n_y, dtype=numpy.float64)
		for j in range(n_y):
			w[j] = Y_counts[j] / ys
	else:
		w = w_pre

	# The (i_min, bin_scale) class of this query in the binned-stage cache.
	h = -1
	if H_keys is not None:
		n_h = H_keys[0, 0]
		for u in range(n_h):
			if H_keys[u+1, 0] == i_min and H_keys[u+1, 1] == bin_scale:
				h = u
				break
		if h == -1 and n_h < H_keys.shape[0] - 1:
			h = n_h
			H_keys[h+1, 0] = i_min
			H_keys[h+1, 1] = bin_scale
			H_keys[0, 0] = n_h + 1

	# Convert the distances to bins and record the histogram of counts. The
	# bin indices are computed in their own loop, which vectorizes, and the
	# scatter-add then runs in the original order. Every index lies in
	# [0, n_bins], because f has n_bins + 1 columns, so int32 holds it exactly.
	#
	# Columns without a cached row, for a 4-letter alphabet, are binned four
	# at a time; consecutive ones have consecutive k, so their four stores
	# land in one gamma_int row, and go through `_binned_block4c`.
	# Of up to three left over, a pair goes through `_binned_block2` and a
	# single one through `_binned_column`.
	zb = numpy.empty(n_y, dtype=numpy.int32)
	zb4 = numpy.empty((4, n_y), dtype=numpy.int32)
	n_todo = 0
	for i in range(nq):
		k = nq - i - 1
		mi = medians[i]
		row = g[i]
		slot = -1
		if q_slot is not None:
			slot = q_slot[i + nq_csum]
			if slot >= 0:
				row = G_cache[slot]

		if slot < 0 and n_a == 4:
			todo[n_todo] = i
			n_todo += 1
			if n_todo == 4:
				i0, i1, i2, i3 = todo[0], todo[1], todo[2], todo[3]
				if i1 == i0 + 1 and i2 == i0 + 2 and i3 == i0 + 3:
					_binned_block4c(g[i0], g[i1], g[i2], g[i3], medians[i0],
						medians[i1], medians[i2], medians[i3], bin_scale,
						offset, w, zb4, gamma_int, uint64(nq - i3 - 1), f[i0],
						f[i1], f[i2], f[i3])
				else:
					_binned_block4(g[i0], g[i1], g[i2], g[i3], medians[i0],
						medians[i1], medians[i2], medians[i3], bin_scale,
						offset, w, zb4, gamma_int, nq - i0 - 1, nq - i1 - 1,
						nq - i2 - 1, nq - i3 - 1, f[i0], f[i1], f[i2], f[i3])
				n_todo = 0
			continue

		hs = -1
		if H_keys is not None:
			if h >= 0 and slot >= 0 and slot < H_int.shape[1]:
				hs = slot
				if H_filled[h, hs]:
					_binned_load(H_int[h, hs], H_f[h, hs], gamma_int, f[i], k)
					continue

		_binned_column(row, mi, bin_scale, offset, w, zb, gamma_int, k, f[i])

		if H_keys is not None and hs >= 0:
			_binned_save(H_int[h, hs], H_f[h, hs], gamma_int, f[i], k)
			H_filled[h, hs] = True

	if n_todo >= 2:
		i0, i1 = todo[0], todo[1]
		_binned_block2(g[i0], g[i1], medians[i0], medians[i1], bin_scale,
			offset, w, zb4, gamma_int, nq - i0 - 1, nq - i1 - 1, f[i0], f[i1])
		todo[0] = todo[2]
		n_todo -= 2
	for u in range(n_todo):
		i = todo[u]
		_binned_column(g[i], medians[i], bin_scale, offset, w, zb, gamma_int,
			nq - i - 1, f[i])

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


@njit(cache=True, inline='always')
def _pairwise_max_packed(x, y, y_csum, z, L, H, x_lo, x_hi, y_lo, y_hi, K, X,
	copy, inplace):
	"""`_pairwise_max_support` with y and y_csum stored as their window only.

	y[i - y_lo] and y_csum[i - y_lo] hold bin i for i in [y_lo, y_hi). Below
	the window y and y_csum are zero; above it y is zero and y_csum is K up
	to bin X and 1 from X on, which is what the dense cumulative sum holds
	there (see `_backgrounds_packed`). Every value read is the one the dense
	row holds, read by the same expression, so z is bitwise the same.
	"""

	if copy:
		for i in range(L, y_lo):
			z[i] = 0.0
		for i in range(y_lo, y_hi):
			z[i] = y[i - y_lo]
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
		z[i] = x[i] * y_csum[i - y_lo] + y[i - y_lo] * x_csum - x[i] * y[i - y_lo]

	if x_hi > y_hi:
		s = max(y_hi, z_lo)
		e = min(x_hi, max(X, s))
		for i in range(s, e):
			z[i] = x[i] * K
		for i in range(e, x_hi):
			z[i] = x[i] * 1.0
	else:
		for i in range(max(x_hi, z_lo), y_hi):
			z[i] = y[i - y_lo] * x_csum

	return z_lo, z_hi


@njit(cache=True, inline='always')
def _pm_packed(x, P, n_packed, a, b, z, L, H, a_lo, a_hi, a_off, a_K, a_X,
	x_lo, x_hi, copy, inplace):
	"""`_pm` on the packed rows: z = max(x, A[a, b]) over [L, H), returning
	the support of z."""

	o, lo, hi = a_off[a, b], a_lo[a, b], a_hi[a, b]
	return _pairwise_max_packed(x, P[o:o+hi-lo],
		P[n_packed+o:n_packed+o+hi-lo], z, L, H, x_lo, x_hi, lo, hi,
		a_K[a, b], a_X[a, b], copy, inplace)


@njit(cache=True, inline='always')
def _pm_fused2_packed(x, y0, c0, p0, q0, K0, X0, z0, y1, c1, p1, q1, K1, X1,
	z1, lo, hi):
	"""`_pm_fused2` with y0, c0, y1 and c1 stored as their windows [p0, q0)
	and [p1, q1) (see `_pairwise_max_packed`).

	[lo, hi) is cut where a window starts or ends or where a constant part
	of a cumulative sum changes. A piece inside both windows runs the loop
	of `_pm_fused2` over views; any other piece runs one loop that uses, for
	each step, its window or its constant. (Specialized loops for the three
	other cases were slower, iteration 100.) Outside its window a step's y is zero and its cumulative sum a constant
	C, so its expression `v * C + y * s - v * y` is `v * C + 0 - 0`, and
	since v is finite and non-negative that is exactly `v * C`. Each running
	sum still takes every v in order, so every value is bitwise the same.
	"""

	s0, s1 = 0.0, 0.0
	i = lo
	while i < hi:
		e = hi
		if i < p0 and p0 < e:
			e = p0
		if i < q0 and q0 < e:
			e = q0
		if i < X0 and X0 < e:
			e = X0
		if i < p1 and p1 < e:
			e = p1
		if i < q1 and q1 < e:
			e = q1
		if i < X1 and X1 < e:
			e = X1

		in0 = i >= p0 and i < q0
		in1 = i >= p1 and i < q1
		C0 = 0.0 if i < p0 else (K0 if i < X0 else 1.0)
		C1 = 0.0 if i < p1 else (K1 if i < X1 else 1.0)

		if in0 and in1:
			xv, u0, u1 = x[i:e], z0[i:e], z1[i:e]
			a0, b0 = y0[i-p0:e-p0], c0[i-p0:e-p0]
			a1, b1 = y1[i-p1:e-p1], c1[i-p1:e-p1]
			for t in range(e - i):
				v = xv[t]
				s0 += v
				t0 = v * b0[t] + a0[t] * s0 - v * a0[t]
				u0[t] = t0
				s1 += t0
				u1[t] = t0 * b1[t] + a1[t] * s1 - t0 * a1[t]
		else:
			for t in range(i, e):
				v = x[t]
				s0 += v
				k0, k1 = uint64(t - p0), uint64(t - p1)
				if in0:
					t0 = v * c0[k0] + y0[k0] * s0 - v * y0[k0]
				else:
					t0 = v * C0
				z0[t] = t0
				s1 += t0
				if in1:
					z1[t] = t0 * c1[k1] + y1[k1] * s1 - t0 * y1[k1]
				else:
					z1[t] = t0 * C1

		i = e


@njit(cache=True, inline='always')
def _pm_chain_packed(B, src, sa, sb, dst, n_steps, P, n_packed, L, H, a_lo,
	a_hi, a_off, a_K, a_X, x_lo, x_hi):
	"""`_pm_chain` on the packed rows, which are never all zero."""

	k, row = 0, src
	while k + 2 <= n_steps:
		lo, hi, z_lo, z_hi = x_lo, x_hi, x_lo, x_hi
		for r in range(k, k+2):
			a, b = sa[r], sb[r]
			hi = max(hi, a_hi[a, b])
			z_lo, z_hi = max(z_lo, a_lo[a, b]), max(z_hi, a_hi[a, b])

		a0, b0, a1, b1 = sa[k], sb[k], sa[k+1], sb[k+1]
		o0, w0 = a_off[a0, b0], a_hi[a0, b0] - a_lo[a0, b0]
		o1, w1 = a_off[a1, b1], a_hi[a1, b1] - a_lo[a1, b1]
		_pm_fused2_packed(B[row], P[o0:o0+w0],
			P[n_packed+o0:n_packed+o0+w0], a_lo[a0, b0], a_hi[a0, b0],
			a_K[a0, b0], a_X[a0, b0], B[dst[k]], P[o1:o1+w1],
			P[n_packed+o1:n_packed+o1+w1], a_lo[a1, b1], a_hi[a1, b1],
			a_K[a1, b1], a_X[a1, b1], B[dst[k+1]], lo, hi)

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
		x_lo, x_hi = _pm_packed(B[row], P, n_packed, sa[k], sb[k], B[dst[k]],
			L, H, a_lo, a_hi, a_off, a_K, a_X, x_lo, x_hi, False,
			dst[k] == row)
		row = dst[k]
		k += 1

	return x_lo, x_hi


@njit(cache=True)
def _backgrounds_packed(f, f_lo, f_hi, B, P, n_packed, nq, n_bins, t_max,
	offset, needed, a_lo, a_hi, L, H, base):
	"""The span backgrounds, their cumulative sums and the rows of B before
	the survival pass, with each span row stored as its nonzero window only.

	Row A[i, j] is nonzero only in [a_lo[i, j], a_hi[i, j]), which the
	caller has checked lies in [0, n) for every row and is never empty. It
	is stored in P[a_off[i, j]:], one window after another, and its
	cumulative sum at the same place plus n_packed. The dense code sums A
	over [c + m, c + m*n_bins], a superset of the window, from 0.0, so its
	cumulative sum is 0 below the window, the window's running sum within
	it, and the total K above it up to bin X = n_bins*(j+1) + c, where it is
	overwritten with 1 from X on. X can be the window's last bin, so the 1
	is written into the window there. The readers use 0, K and 1 outside the
	window. The convolution adds the same terms in the same order into a
	row cleared over its window, which is exactly where it writes, so every
	stored value is bitwise the dense one.

	Column c of a row of B holds bin c + base. The B build reads and writes
	only bins in [L, H), and depends on bins only through their differences
	and their order relative to the window ends and X, so it runs unchanged
	on bins shifted down by base.
	"""

	a_off = numpy.empty((nq, nq), dtype='int64')
	a_K = numpy.empty((nq, nq), dtype='float64')
	a_X = numpy.empty((nq, nq), dtype='int64')
	p = int64(0)
	for i in range(nq):
		for j in range(i, nq):
			a_off[i, j] = p
			p += a_hi[i, j] - a_lo[i, j]
			a_X[i, j] = int64(n_bins) * int64(j+1) + int64(offset) * int64(nq
				- j + i - 1)

	for i in range(nq):
		for j in range(i, nq):
			o, w = a_off[i, j], a_hi[i, j] - a_lo[i, j]
			row = P[o:o+w]
			if i == j:
				for l in range(f_lo[j], f_hi[j]+1):
					row[l - f_lo[j]] = f[j, l]
			else:
				row[:] = 0

				# In the dense row, bin k + c + offset of A[i, j-1] times
				# f[j, l_lo + s] is added into bin k + c + l_lo + s of A[i, j],
				# for k over the support [k_lo, k_hi] of A[i, j-1] in the
				# `k` coordinate. Here both are shifted by their windows.
				c = int64(offset) * int64(nq - j + i - 1)
				k_lo = a_lo[i, j-1] - c - int64(offset)
				k_hi = a_hi[i, j-1] - 1 - c - int64(offset)
				po = a_off[i, j-1]
				prev = P[po:po + a_hi[i, j-1] - a_lo[i, j-1]]
				l_lo, l_hi = f_lo[j], f_hi[j]
				_convolve_span(prev, row, f[j, l_lo:l_hi+1], k_lo, k_hi,
					-k_lo, c + l_lo - a_lo[i, j])

	for i in range(nq):
		for j in range(i, nq):
			o, w = a_off[i, j], a_hi[i, j] - a_lo[i, j]
			a, s = P[o:o+w], P[n_packed+o:n_packed+o+w]
			acc = 0.0
			for k in range(w):
				acc += a[k]
				s[k] = acc

			a_K[i, j] = acc
			for k in range(max(a_X[i, j] - a_lo[i, j], 0), w):
				s[k] = 1.0

	# From here bins are B's columns: bin - base.
	c_lo = numpy.empty((nq, nq), dtype='int64')
	c_hi = numpy.empty((nq, nq), dtype='int64')
	for i in range(nq):
		for j in range(i, nq):
			c_lo[i, j] = a_lo[i, j] - base
			c_hi[i, j] = a_hi[i, j] - base
			a_X[i, j] -= base
	Lc, Hc = L - base, H - base

	# The B build of `_backgrounds_dense`, reading the packed rows. Pass 0
	# builds the chain for rows nq..t_max and pass r > 0 row r; both end in
	# the one call to `_pm_chain_packed`, which is inlined, so it is 
	# compiled once.
	n_steps_max = 2*nq + t_max
	sa = numpy.empty(n_steps_max, dtype='int64')
	sb = numpy.empty(n_steps_max, dtype='int64')
	dst = numpy.empty(n_steps_max, dtype='int64')

	for r in range(min(nq, t_max+1)):
		ns = 0
		if r == 0:
			if t_max < nq:
				continue

			b_lo, b_hi, src = Lc, Hc, 1
			if nq > 1:
				b_lo, b_hi = _pm_packed(B[0], P, n_packed, 0, 0, B[1], Lc, Hc, 
					c_lo, c_hi, a_off, a_K, a_X, b_lo, b_hi, True, False)
				sa[ns], sb[ns], dst[ns] = nq-1, nq-1, 1
				ns += 1

			for i in range(2, nq):
				sa[ns], sb[ns], dst[ns] = 0, i-1, i
				sa[ns+1], sb[ns+1], dst[ns+1] = nq-i, nq-1, i
				ns += 2

			for i in range(nq, t_max+1):
				if i == 1:
					b_lo, b_hi = _pm_packed(B[0], P, n_packed, 0, nq-1, B[1], 
						Lc, Hc, c_lo, c_hi, a_off, a_K, a_X, b_lo, b_hi, True, 
						False)
					continue
				sa[ns], sb[ns], dst[ns] = 0, nq-1, i
				ns += 1
		else:
			if needed is not None and not needed[r]:
				continue

			b_lo, b_hi, src = Lc, Hc, r
			b_lo, b_hi = _pm_packed(B[r], P, n_packed, 0, r-1, B[r], Lc, Hc, 
				c_lo, c_hi, a_off, a_K, a_X, Lc, Hc, True, True)

			for j in range(1, nq - r + 1):
				sa[ns], sb[ns], dst[ns] = j, j+r-1, r
				ns += 1

			for j in range(r-1):
				sa[ns], sb[ns], dst[ns] = 0, j, r
				sa[ns+1], sb[ns+1], dst[ns+1] = nq-1-j, nq-1, r
				ns += 2

		b_lo, b_hi = _pm_chain_packed(B, src, sa, sb, dst, ns, P, n_packed, Lc,
			Hc, c_lo, c_hi, a_off, a_K, a_X, b_lo, b_hi)


@njit(cache=True)
def _backgrounds_dense(f, f_lo, f_hi, A, B, A_csum, nq, n_bins, t_max, offset,
	needed, a_lo, a_hi, L, H, n):
	"""The span backgrounds A, their cumulative sums A_csum and the rows of B
	before the survival pass, in the dense layout: every row of A and A_csum
	is n bins long. `_p_value_backgrounds` uses this when a query has an
	empty column or a support that reaches n."""

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
				for l in range(n_bins+1):
					l = uint64(l)
					A[i, j, l+c] = f[j, l]

				k_lo, k_hi = f_lo[j], f_hi[j]
			else:
				l_lo, l_hi = f_lo[j], f_hi[j]

				_convolve_span(A[i, j-1], A[i, j], f[j, l_lo:l_hi+1], max(k_lo, 0),
					min(k_hi, numpy.int64(n_bins*j)), int64(c+offset), int64(c)+l_lo)

				k_lo, k_hi = k_lo + l_lo, k_hi + l_hi


	_A_cumsum(A, A_csum, nq, n_bins, offset, n)

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



@njit(cache=True)
def _p_value_backgrounds_windowed(f, A, Bf, A_csum, nq, n_bins, t_max, offset,
	needed, packed):
	"""An internal function that calculates the backgrounds for p-values.

	This method takes in the histogram of integerized scores `f` and returns 
	the background probabilities of each overlap achieving a given score. 
	These scores are calculated for the complete overlap of the query and
	target, but also for all overhangs where only part of the query and the
	target are overlapping (on either end). Additionally, background
	probabilities are calculated for all spans across the query for when the
	target is smaller than the query and has to be scanned against it.

	Only the rows `t` with `needed[t]` true are finished; the others hold 
	unspecified values. `needed` may be None, meaning every row.

	The rows of B are stored in the flat array `Bf` as t_max+1 rows of S 
	entries, and the function returns (lo, S). The survival value of a
	target of length t at score s > 0, which the dense layout holds at 
	B[t, s-1], is Bf[t*S + k] with k = s - lo clipped to [0, S-1] 
	(`_b_index`). Normally a row holds only the bins [L, H) where the
	survival function varies, with the constant 1.0 of the bins below L in
	front and the constant of the bins from H on behind: S = H - L + 2 and
	lo = L. A query that takes the dense path below, and any call that
	needs row 0, which is not constant outside [L, H), gets the dense 
	layout: S = n and lo = 1, so k = s - 1.

	The span backgrounds A[i, j] and their cumulative sums are stored packed,
	each row as its nonzero window only (see `_backgrounds_packed`), in
	`packed` when it is given, and A and A_csum are then not written.
	Without `packed` they are stored in the memory of A_csum, and each A[i, i]
	is set to the distribution of column i. A query with an empty column, a
	support that reaches n, or too little room uses the dense layout in A and
	A_csum (`_backgrounds_dense`).
	"""

	n = n_bins*nq + nq*offset

	# First and last nonzero bin of each query column's histogram. Every term
	# is non-negative, so a skipped zero term would only have added +0.0, and
	# the remaining terms are still added in the same order: bitwise-exact.
	f_lo = numpy.empty(nq, dtype='int64')
	f_hi = numpy.empty(nq, dtype='int64')
	for j in range(nq):
		f_lo[j], f_hi[j] = n_bins+1, 0
		for l in range(n_bins+1):
			if f[j, l] != 0:
				f_hi[j] = l
				if f_lo[j] > n_bins:
					f_lo[j] = l

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

	# The packed layout needs every support inside [0, n) and no empty column.
	dense = H > int64(n)
	for j in range(nq):
		if f_lo[j] > n_bins:
			dense = True

	if H > int64(n):
		H = int64(n)
	if H <= L:
		L, H, dense = 0, 0, True

	# A[a, b] is zero outside [a_lo, a_hi) (clipped to [L, H)), and a_lo is 
	# -1 when the row is all zero; its A_csum is then not zero below the fill
	# of 1, so it goes through the unrestricted `_pairwise_max_window` and 
	# the result's support is taken to be all of [L, H).
	for i in range(nq):
		for j in range(i, nq):
			if a_lo[i, j] >= 0:
				a_hi[i, j] = min(a_hi[i, j], H)

	# Row 0 is not constant outside [L, H), so only the dense layout holds it.
	if needed is None or needed[0]:
		dense = True

	lo, n_col = int64(1), int64(n)
	n_packed = int64(0)
	if packed is None:
		P = A_csum.reshape(A_csum.size)
	else:
		P = packed

	if not dense:
		for i in range(nq):
			for j in range(i, nq):
				n_packed += a_hi[i, j] - a_lo[i, j]

		dense = 2*n_packed > P.shape[0] or \
			(int64(t_max) + 1) * (H - L + 2) > Bf.shape[0]

	if dense:
		B = Bf[:(t_max+1)*n].reshape((t_max+1, n))

		# Never taken on JASPAR or on random motifs of widths 1-40, so
		# `_backgrounds_dense` is called through object mode and compiled on
		# first use rather than with every caller (3.5 s of a cold first
		# call, plus its code relinked into each enclosing library). The
		# scalars travel as one-element arrays so the call is typed with
		# their numba types (offset is uint64 from `_tomtom`).
		s_nq, s_nb, s_tm = numpy.full(1, nq), numpy.full(1, n_bins), \
			numpy.full(1, t_max)
		s_off, s_L, s_H, s_n = numpy.full(1, offset), numpy.full(1, L), \
			numpy.full(1, H), numpy.full(1, n)
		with numba.objmode():
			_backgrounds_dense(f, f_lo, f_hi, A, B, A_csum, s_nq[0], s_nb[0],
				s_tm[0], s_off[0], needed, a_lo, a_hi, s_L[0], s_H[0], s_n[0])

		# Row 0 is the all -1 starting point and is never a real distribution.
		if needed is None or needed[0]:
			for j in range(n):
				B[0, j] = -1
			for j in range(1, n):
				B[0, j] += B[0, j-1]
			for j in range(n):
				b = 1 - B[0, j]
				B[0, j] = b if b > 0 else 0.0
	else:
		if packed is None:
			c = int64(offset) * int64(nq - 1)
			for i in range(nq):
				A[i, i] = 0
				for l in range(n_bins+1):
					A[i, i, l+c] = f[i, l]

		# Column c of a row holds bin c + L - 1: the survival pass below runs
		# over columns [1, H - L + 1), and the constants go in the two ends.
		n_col = H - L + 2
		B = Bf[:(t_max+1)*n_col].reshape((t_max+1, n_col))
		_backgrounds_packed(f, f_lo, f_hi, B, P, n_packed, nq, n_bins, t_max,
			offset, needed, a_lo, a_hi, L, H, L - 1)

		lo, L, H = L, int64(1), n_col - 1

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
		for j in range(H, n_col):
			B[i, j] = tail

	return lo, n_col


@njit(cache=True, inline='always')
def _b_index(score, lo, n_col):
	"""The column of a row of B holding the survival value at `score` > 0,
	for a layout (lo, n_col) returned by `_p_value_backgrounds_windowed`."""

	return uint64(min(max(int64(score) - lo, int64(0)), n_col - 1))


@njit(cache=True)
def _p_value_backgrounds(f, A, B, A_csum, nq, n_bins, t_max, offset, 
	needed=None, packed=None):
	"""`_p_value_backgrounds_windowed` with B returned in the dense layout:
	B[t, j] is the survival value at score j + 1 for j < n, for every 
	finished row t, and row 0 is finished when `needed` is None or 
	needed[0]. B holds t_max+1 rows of at least n entries.
	"""

	n = n_bins*nq + nq*offset

	nd = numpy.ones(t_max+1, dtype=numpy.bool_)
	if needed is not None:
		nd[:] = needed[:t_max+1]
	row0 = nd[0]
	nd[0] = False

	Bf = numpy.empty((t_max+1) * (n+2), dtype=B.dtype)
	lo, n_col = _p_value_backgrounds_windowed(f, A, Bf, A_csum, nq, n_bins,
		t_max, offset, nd, packed)

	for i in range(1, t_max+1):
		if nd[i]:
			for j in range(n):
				B[i, j] = Bf[uint64(i*n_col) + _b_index(j+1, lo, n_col)]

	if row0:
		for j in range(n):
			B[0, j] = -1
		for j in range(1, n):
			B[0, j] += B[0, j-1]
		for j in range(n):
			b = 1 - B[0, j]
			B[0, j] = b if b > 0 else 0.0
			

@njit(cache=True)
def _p_values(gamma, B_cdfs, rr_inv, T_lens, iq, nq, offset, results,
	reverse_complement=1, b_lo=1):
	"""An internal function for calculating the best match and p-values.

	Chooses the width of the running sums and calls `_p_values_sums`. Each
	sum is nq * offset plus at most nq int16 values of `gamma`, so it lies in
	[-nq * 32768, nq * (offset + 32767)]. When that fits in int32 the sums
	are exact in int32, which halves the accumulation's vector width;
	otherwise they are kept in int64.

	t_sums has 8 slots past the longest scan and is filled with the dtype's
	minimum, which no sum can reach; `_p_values_sums` reads it back from the
	last slot as its padding value.

	The int64 sums need nq * (offset + 32768) > 2**31, a query of more than
	65,535 columns, so they run in object mode (`_p_values_int64`) and are
	compiled on first use rather than with every caller.
	"""

	# Sized by the longest target, not by gamma, whose rows are the unique
	# target columns and can be fewer than a target's length after hashing.
	n_sums = T_lens.max() + nq - 1 + 8

	if int64(nq) * (int64(offset) + 32768) <= 2147483647:
		t_sums = numpy.full(n_sums, -2147483648, dtype='int32')
		_p_values_sums(gamma, B_cdfs, rr_inv, T_lens, iq, nq, offset, 
			results, reverse_complement, t_sums, b_lo)
	else:
		# The scalars travel as one-element arrays so the object-mode call
		# is typed with their numba types (offset is uint64 from `_tomtom`).
		# `_p_values_sums` only tests reverse_complement == 1.
		s_iq, s_nq, s_off = numpy.full(1, iq), numpy.full(1, nq), \
			numpy.full(1, offset)
		s_rc = numpy.full(1, int64(1 if reverse_complement == 1 else 0))
		s_lo = numpy.full(1, b_lo)
		with numba.objmode():
			_p_values_int64(gamma, B_cdfs, rr_inv, T_lens, s_iq[0], s_nq[0],
				s_off[0], results, s_rc[0], n_sums, s_lo[0])


def _p_values_int64(gamma, B_cdfs, rr_inv, T_lens, iq, nq, offset, results,
	reverse_complement, n_sums, b_lo):
	"""The int64 branch of `_p_values`, called from its object-mode block."""

	t_sums64 = numpy.full(n_sums, -9223372036854775807 - 1, dtype='int64')
	_p_values_sums(gamma, B_cdfs, rr_inv, T_lens, iq, nq, offset, results,
		reverse_complement, t_sums64, b_lo)


@njit(cache=True)
def _p_values_sums(gamma, B_cdfs, rr_inv, T_lens, iq, nq, offset, results,
	reverse_complement, t_sums, b_lo):
	"""An internal function for calculating the best match and p-values.

	This function will take in the integerized score matrix `gamma` and
	background distributions `B_cdfs` and calculate the best overlap.
	The best overlap is calculated as the best sum of scores across the
	alignment, minus a penalty for each unaligned column. After finding
	a new best overlap, the p-value is calculated by comparing the
	score to the background distribution.

	Targets 0..iq are skipped, and so are their reverse complements, which
	start at len(T_lens) // 2 only when `reverse_complement` is 1.

	The p-value at score s > 0 is B_cdfs[nt, _b_index(s, b_lo, n_col)] for
	B_cdfs's rows of n_col entries (see `_p_value_backgrounds_windowed`);
	b_lo = 1 is the dense layout, B_cdfs[nt, s-1].
	"""

	b_lo = int64(b_lo)
	n_col = int64(B_cdfs.shape[1])

	n = len(T_lens) // 2 if reverse_complement == 1 else len(T_lens)
	total_offset = uint64(0)

	# The dtype's minimum, which `_p_values` put in the last slot and which
	# no sum can equal.
	pad = t_sums[len(t_sums)-1]

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
		# The sums past m are set to `pad`, so both passes below run over a
		# multiple of 8 and, for int32 sums, take the 8-lane path with no
		# scalar tail.
		# `pad` is below every sum, so it changes neither M nor kf/kl.
		m = int64(nt) + int64(nq) - 1
		for j in range(8):
			t_sums[m+j] = pad
		mr = (m + 7) // 8 * 8

		M = pad
		for k in range(mr):
			M = max(M, t_sums[k])

		# The first and last positions holding M, by branchless integer
		# min/max, which vectorize; every position equal to M lies in
		# kf..kl and is still visited in increasing k.
		kf32 = int32(mr)
		kl32 = int32(0)
		for k in range(mr):
			e = t_sums[k] == M
			kf32 = min(kf32, int32(k) if e else int32(mr))
			kl32 = max(kl32, int32(k) if e else int32(0))
		kf = int64(kf32)
		kl = int64(kl32)

		# One position holds a positive M: the visit below would pass the
		# first test (results[i, 1] is 0) and not the tie test, and write
		# all four fields from it.
		if kf == kl and M > 0:
			results[i, 0] = B_cdfs[nt, _b_index(M, b_lo, n_col)]
			results[i, 1] = M
			results[i, 2] = kf - nq + 1
			results[i, 3] = min(kf+1, nq) - max(0, kf-int64(nt)+1)
			total_offset += nt
			continue

		for k in range(kf, kl+1):
			score = t_sums[k]
			if score != M:
				continue

			overlap = min(k+1, nq) - max(0, k-nt+1)
			if score >= results[i, 1]:
				if score == results[i, 1] and results[i, 2] >= overlap:
					continue

				results[i, 0] = B_cdfs[nt, _b_index(score, b_lo, n_col)] if score > 0 else 1.0
				results[i, 1] = score
				results[i, 2] = k - nq + 1
				results[i, 3] = overlap

		total_offset += nt


@njit(cache=True)
def _p_values_batch(G, Bs, rr_inv, T_lens, nq, offsets, active, results, nb,
	s_proto, b_lo, b_col):
	"""Allocates the running sums and calls `_p_values_batch_sums`.

	A full batch passes NB = 16 as a constant, so numba compiles a kernel
	specialized to it; batches of 2, 4 and 8 share one kernel with NB a
	runtime value. A kernel per batch size tied on runtime for 16 only and
	cost about 2 s of cold compile (iteration 113); a runtime NB for all
	sizes was 0.01 s slower.

	The sums, the fill and the tracked maxima take the dtype of `s_proto`:
	int16 for an int8 `G`, int32 for an int16 `G`. `_tomtom` puts a query in
	the batch only when its partial sums fit that dtype (`_sums_fit`).
	Inactive queries get a fill of 0, since their lanes are never read.
	"""

	n_fill = (uint64(T_lens.max()) + uint64(nq) - 1) * uint64(nb)
	t_sums = numpy.empty(n_fill + 16, dtype=s_proto.dtype)
	fill = numpy.empty(n_fill, dtype=s_proto.dtype)
	mv = numpy.empty(16, dtype=s_proto.dtype)
	m_init = -(int64(1) << int64(8 * s_proto.itemsize - 1))
	if nb == 16:
		_p_values_batch_sums(G, Bs, rr_inv, T_lens, nq, offsets, active,
			results, t_sums, fill, mv, m_init, 16, b_lo, b_col)
	else:
		_p_values_batch_sums(G, Bs, rr_inv, T_lens, nq, offsets, active,
			results, t_sums, fill, mv, m_init, nb, b_lo, b_col)


@njit(cache=True, inline='always')
def _sums_fit(nq, offset, g_max, s_max):
	"""Whether a query's partial sums fit in [-s_max-1, s_max]: they start
	at nq * offset and receive at most nq values of `gamma_int`, each in
	[-g_max-1, g_max], whatever its contents."""

	lo = int64(nq) * (int64(offset) - int64(g_max) - 1)
	hi = int64(nq) * (int64(offset) + int64(g_max))
	return lo >= -int64(s_max) - 1 and hi <= int64(s_max)


@njit(cache=True)
def _p_values_batch_sums(G, Bs, rr_inv, T_lens, nq, offsets, active, results,
	t_sums, fill, mv, m_init, NB, b_lo, b_col):
	"""`_p_values` for NB queries of the same width nq, in one target loop.

	NB is 16 or the batch size (see `_p_values_batch`). Row r of `G` holds
	the NB queries' `gamma_int` rows side by side, as
	G[r, l*NB + b] = gamma_int_b[r, l], and the running
	sums are interleaved the same way: position p of query b is
	t_sums[p*NB + b]. Row k of a target then adds one contiguous run of
	nq*NB values into t_sums[k*NB : (k+nq)*NB], so the gather through
	`rr_inv` and the loop setup are paid once for NB queries. Each query's
	sums receive the same integer terms as in `_p_values_sums`, and every
	active query's sums must fit the dtype of `t_sums` (`_sums_fit`).

	Each query's maximum M and the first and last positions holding it are
	tracked over the target's positions once its rows are added: the same M,
	kf and kl as the separate scans in `_p_values_sums`.
	That update always runs over 16 lanes, which LLVM vectorizes with
	selects; over NB < 16 lanes the selects became data-dependent branches.
	Lanes at or above NB read the next position and are ignored, so t_sums
	has 16 spare entries. The tie-rule scan over kf..kl, the p-value from
	Bs[b] and the outputs are the per-query ones, in the same order.
	Queries with active[b] False are accumulated but not scanned, and their
	`results[b]` is left alone. There is no target skipping (iq is -1).

	`t_sums` (n_fill + 16 entries), `fill` (n_fill) and `mv` (16) share one
	integer dtype, chosen by `_tomtom`; `m_init` is its minimum.

	Bs[b] holds query b's rows of B flat, b_col[b] entries each, in the
	layout (b_lo[b], b_col[b]) that `_p_value_backgrounds_windowed` returned.
	The lookup clips the column only from below: a score is at most 
	nq * n_score_bins <= n, the last column of the dense layout, and in the
	windowed layout, which a query gets only when no column's histogram is
	empty, at most H - 1, the largest score any overlap span's background
	supports (every value in gamma_int has nonzero weight in f).
	"""

	nb = uint64(NB)
	w = uint64(nq) * nb
	n_fill = uint64(len(fill))

	# Allocated here and in `_p_values_batch`, not as per-thread rows of a
	# shared array: mv, fv and lv are written at every position, and
	# neighbouring threads' rows sharing a cache line made 8 threads up to
	# 25x slower.
	fv = numpy.empty(16, dtype=numpy.int32)
	lv = numpy.empty(16, dtype=numpy.int32)
	for p in range(n_fill // nb):
		for b in range(NB):
			if active[b]:
				fill[uint64(p)*nb + uint64(b)] = nq * offsets[b]
			else:
				fill[uint64(p)*nb + uint64(b)] = 0

	total_offset = uint64(0)
	for i, nt in enumerate(T_lens):
		nt = uint64(nt)
		m = nt + uint64(nq) - 1
		n = m * nb
		for k in range(n):
			k = uint64(k)
			t_sums[k] = fill[k]

		for b in range(16):
			mv[b] = m_init
			fv[b] = 0
			lv[b] = 0

		for k in range(m):
			k = uint64(k)
			base = k * nb
			if k < nt:
				r = uint64(rr_inv[total_offset + k])
				for j in range(w):
					j = uint64(j)
					t_sums[base + j] += G[r, j]

		# The tracking runs after the target's rows are added, four positions
		# per pass over the 16 lanes, so mv, fv and lv are loaded and stored
		# once per four positions instead of once per position; each lane
		# still sees the positions in increasing k with the same comparisons.
		k = uint64(0)
		while k + uint64(4) <= m:
			base = k * nb
			k1, k2, k3 = k + uint64(1), k + uint64(2), k + uint64(3)
			for b in range(16):
				ub = uint64(b)
				mb, fb, lb = mv[b], fv[b], lv[b]
				v = t_sums[base + ub]
				fb = numpy.int32(k) if v > mb else fb
				lb = numpy.int32(k) if v >= mb else lb
				mb = max(mb, v)
				v = t_sums[base + nb + ub]
				fb = numpy.int32(k1) if v > mb else fb
				lb = numpy.int32(k1) if v >= mb else lb
				mb = max(mb, v)
				v = t_sums[base + uint64(2) * nb + ub]
				fb = numpy.int32(k2) if v > mb else fb
				lb = numpy.int32(k2) if v >= mb else lb
				mb = max(mb, v)
				v = t_sums[base + uint64(3) * nb + ub]
				fb = numpy.int32(k3) if v > mb else fb
				lb = numpy.int32(k3) if v >= mb else lb
				mb = max(mb, v)
				mv[b], fv[b], lv[b] = mb, fb, lb
			k += uint64(4)

		while k < m:
			base = k * nb
			for b in range(16):
				v = t_sums[base + uint64(b)]
				mb = mv[b]
				fv[b] = numpy.int32(k) if v > mb else fv[b]
				lv[b] = numpy.int32(k) if v >= mb else lv[b]
				mv[b] = max(mb, v)
			k += uint64(1)

		for b in range(NB):
			if not active[b]:
				continue

			res = results[b]
			B_cdfs = Bs[b]
			res[i, 0] = 1
			res[i, 1] = 0
			lo_b, col_b = b_lo[b], b_col[b]
			row_b = uint64(nt) * uint64(col_b)

			ub = uint64(b)
			M = mv[b]
			kf = int64(fv[b])
			kl = int64(lv[b])
			for k in range(kf, kl+1):
				score = t_sums[uint64(k)*nb + ub]
				if score != M:
					continue

				overlap = min(k+1, nq) - max(0, k-nt+1)
				if score >= res[i, 1]:
					if score == res[i, 1] and res[i, 2] >= overlap:
						continue

					res[i, 0] = B_cdfs[row_b + uint64(max(int64(score) - lo_b, 
						int64(0)))] if score > 0 else 1.0
					res[i, 1] = score
					res[i, 2] = k - nq + 1
					res[i, 3] = overlap

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


@njit(cache=True)
def _integer_histogram_int16(Y, gamma, gamma_int, f, medians, Y_counts,
	nq_csum, nq, i_min, bin_scale, offset, q_slot, G_cache, H_keys, H_filled,
	H_int, H_f, w_pre):
	"""`_integer_histogram` into an int16 `gamma_int`, called through object
	mode so that it compiles on first use. `_tomtom` needs it only for a
	query whose offset does not fit the int8 batch. The scalars travel as
	one-element arrays so the call is typed with their numba types."""

	s_c, s_nq, s_im = numpy.full(1, nq_csum), numpy.full(1, nq), \
		numpy.full(1, i_min)
	s_bs, s_off = numpy.full(1, bin_scale), numpy.full(1, offset)
	with numba.objmode():
		_integer_histogram(Y, gamma, gamma_int, f, medians, Y_counts, s_c[0],
			s_nq[0], s_im[0], s_bs[0], s_off[0], q_slot, G_cache, H_keys,
			H_filled, H_int, H_f, w_pre)


@njit(cache=True, inline='always')
def _transpose8(r0, r1, r2, r3, r4, r5, r6, r7):
	"""Transpose the 8x8 byte matrix whose row b is word r_b (byte p at
	bits 8p): word p of the result holds byte p of every r_b, r_b's at bits
	8b. Three rounds of masked swaps."""

	m32 = uint64(0x00000000FFFFFFFF)
	t = ((r0 >> uint64(32)) ^ r4) & m32
	r0 ^= t << uint64(32)
	r4 ^= t
	t = ((r1 >> uint64(32)) ^ r5) & m32
	r1 ^= t << uint64(32)
	r5 ^= t
	t = ((r2 >> uint64(32)) ^ r6) & m32
	r2 ^= t << uint64(32)
	r6 ^= t
	t = ((r3 >> uint64(32)) ^ r7) & m32
	r3 ^= t << uint64(32)
	r7 ^= t

	m16 = uint64(0x0000FFFF0000FFFF)
	t = ((r0 >> uint64(16)) ^ r2) & m16
	r0 ^= t << uint64(16)
	r2 ^= t
	t = ((r1 >> uint64(16)) ^ r3) & m16
	r1 ^= t << uint64(16)
	r3 ^= t
	t = ((r4 >> uint64(16)) ^ r6) & m16
	r4 ^= t << uint64(16)
	r6 ^= t
	t = ((r5 >> uint64(16)) ^ r7) & m16
	r5 ^= t << uint64(16)
	r7 ^= t

	m8 = uint64(0x00FF00FF00FF00FF)
	t = ((r0 >> uint64(8)) ^ r1) & m8
	r0 ^= t << uint64(8)
	r1 ^= t
	t = ((r2 >> uint64(8)) ^ r3) & m8
	r2 ^= t << uint64(8)
	r3 ^= t
	t = ((r4 >> uint64(8)) ^ r5) & m8
	r4 ^= t << uint64(8)
	r5 ^= t
	t = ((r6 >> uint64(8)) ^ r7) & m8
	r6 ^= t << uint64(8)
	r7 ^= t
	return r0, r1, r2, r3, r4, r5, r6, r7


@njit(cache=True)
def _interleave(Gs, G, n, nb):
	"""Write the nb packed matrices Gs[b, :n] into a batch's interleaved
	matrix, G[p*nb + b] = Gs[b, p] (see `_p_values_batch`).

	`_integer_distances_and_histogram` writes a query's `gamma_int` a few
	columns per sweep over the targets. Into its lane of G directly, every
	sweep touches every cache line of the batch's n*nb matrix, several MB;
	into its own packed row of Gs the writes stay in L2. For one-byte
	values and nb of 8 or 16 the copy moves eight positions of eight lanes
	as eight 64-bit words through a byte transpose, with a scalar loop for
	the rest. Gs's rows must start on 8-byte boundaries, and the words are
	little-endian. Copied byte by byte, the copy cost as much as the L2
	saved (iteration 109). The values are copied, not changed.
	"""

	p0 = uint64(0)
	if G.itemsize == 1 and nb == 16:
		n8 = uint64(n) // uint64(8)
		G64 = G.view(numpy.uint64)
		S = Gs[:16].view(numpy.uint64)
		for q in range(n8):
			q = uint64(q)
			o = q * uint64(16)
			c0, c1, c2, c3, c4, c5, c6, c7 = _transpose8(S[0, q], S[1, q],
				S[2, q], S[3, q], S[4, q], S[5, q], S[6, q], S[7, q])
			d0, d1, d2, d3, d4, d5, d6, d7 = _transpose8(S[8, q], S[9, q],
				S[10, q], S[11, q], S[12, q], S[13, q], S[14, q], S[15, q])
			G64[o], G64[o+uint64(1)] = c0, d0
			G64[o+uint64(2)], G64[o+uint64(3)] = c1, d1
			G64[o+uint64(4)], G64[o+uint64(5)] = c2, d2
			G64[o+uint64(6)], G64[o+uint64(7)] = c3, d3
			G64[o+uint64(8)], G64[o+uint64(9)] = c4, d4
			G64[o+uint64(10)], G64[o+uint64(11)] = c5, d5
			G64[o+uint64(12)], G64[o+uint64(13)] = c6, d6
			G64[o+uint64(14)], G64[o+uint64(15)] = c7, d7
		p0 = n8 * uint64(8)
	elif G.itemsize == 1 and nb == 8:
		n8 = uint64(n) // uint64(8)
		G64 = G.view(numpy.uint64)
		S = Gs[:8].view(numpy.uint64)
		for q in range(n8):
			q = uint64(q)
			o = q * uint64(8)
			c0, c1, c2, c3, c4, c5, c6, c7 = _transpose8(S[0, q], S[1, q],
				S[2, q], S[3, q], S[4, q], S[5, q], S[6, q], S[7, q])
			G64[o], G64[o+uint64(1)], G64[o+uint64(2)] = c0, c1, c2
			G64[o+uint64(3)], G64[o+uint64(4)], G64[o+uint64(5)] = c3, c4, c5
			G64[o+uint64(6)], G64[o+uint64(7)] = c6, c7
		p0 = n8 * uint64(8)

	for p in range(p0, uint64(n)):
		p = uint64(p)
		for b in range(nb):
			G[p*uint64(nb) + uint64(b)] = Gs[b, p]


# Only the explicit prange loop is parallelized. With parallel=True numba
# also turns the whole-array numpy calls outside it (zeros, cumsum, sum) into
# parallel loops, each compiled and linked separately, for arrays of a few
# thousand elements. A ParallelOptions object, not a dict: numba pops the
# dict's keys on the first compile, so a second signature would get every
# option back.
_PRANGE_ONLY = ParallelOptions({'comprehension': False, 'reduction': False,
	'inplace_binop': False, 'setitem': False, 'numpy': False, 
	'stencil': False, 'fusion': False, 'prange': True})


@njit(parallel=_PRANGE_ONLY, cache=True)
def _tomtom(Q, T, Q_lens, T_lens, Q_norm, T_norm, rr_inv, rr_counts, n_nearest, 
	n_score_bins, n_median_bins, n_cache, n_threads, reverse_complement,
	q_slot, q_cached, results, _A, _A_csum, _B, _G, _Gs, direct, G_cache, 
	g_max, s_proto, order):
	"""An internal function implementing the TOMTOM algorithm.

	This internal function is necessary to handle the numba component of the
	implementation. Here, scratchboard memory is allocated for each thread and
	the main parallel loop is called. Additionally, if reverse complements are
	being considered, values are merged across both strands.

	`_G`'s dtype and the empty array `s_proto`'s are the batched
	`gamma_int`'s and its sums': int8 and int16 with `g_max` 127, or int16
	and int32 with `g_max` -1 (see `tomtom`). Passing the dtypes in compiles
	only the batched kernels that the call uses.
	"""

	T_max = max(T_lens)

	# `_p_values` reads B only at the rows that are target lengths.
	needed = numpy.zeros(T_max+1, dtype=numpy.bool_)
	for t in T_lens:
		needed[t] = True
	
	Q_offsets = numpy.zeros(len(Q_lens)+1, dtype='int64')
	for q in range(len(Q_lens)):
		Q_offsets[q+1] = Q_offsets[q] + Q_lens[q]
	Q_max = max(Q_lens)
	
	n_in_targets = len(T_lens) // 2 if reverse_complement else len(T_lens)
	nt = T.shape[-1]

	# Re-usable workspace for each thread instead of re-allocating
	# and freeing large arrays for each example.
	n_len = Q_max*n_score_bins + Q_max*n_cache
	
	_gamma = numpy.empty((n_threads, nt, Q_max), dtype='float64')
	# Flat per thread; each query takes a (nt, nq) view, so a row holds only
	# the nq columns that are written and read, and the rows are packed.
	_gamma_int = numpy.empty((n_threads, nt*Q_max), dtype='int16')

	# Up to n_batch queries of the same width share one `_p_values_batch`
	# target loop. Their `gamma_int` matrices are interleaved in _G (see
	# `_p_values_batch`), and each keeps its own backgrounds and results.
	# A batch is 2, 4, 8 or 16 queries; 16 has its own compiled kernel and
	# the others share one. 16 was faster than 8 or 32 (iteration 74).
	#
	# _G is int8 when g_max is 127: every value x - offset (x in
	# [0, n_score_bins]) then fits once `_fits_gamma` holds, and the batch's
	# sums can be int16 (see `_p_values_batch`). A query whose offset does
	# not fit is rerun into its int16 `_gamma_int` and runs alone.
	n_batch = 16
	_f = numpy.empty((n_threads, Q_max, n_score_bins+1), dtype='float64')

	# A and A_csum are flat per thread; each query takes a contiguous
	# (nq, nq, n) view so its working set is not strided by n_len.
	# `_A`, `_A_csum`, `_B`, `_G`, `G_cache` and `results` are allocated by
	# `tomtom` with numpy; see the comment there.

	_medians = numpy.empty((n_threads, Q_max), dtype='float64')
	_median_bins = numpy.empty((n_threads, n_median_bins, 2), dtype='float64')

	_results = numpy.empty((n_threads, n_batch, len(T_lens), 5), 
		dtype='float64')

	# Query columns that occur more than once have their distance row, min,
	# max and median computed once here; every query reads them from these.
	S_cache = numpy.empty((len(q_cached), 3), dtype='float64')
	_fill_column_cache(Q, T, Q_norm, T_norm, rr_counts, q_cached, G_cache,
		S_cache, n_median_bins)

	# The target weights and the median's half-count depend only on the
	# targets, so they are computed once here, by the expressions
	# `_integer_distances_and_histogram` would use, and read by every query.
	halfway = _halfway(rr_counts)
	ys = numpy.sum(rr_counts)
	# The block medians count in int32 when every count, and so every bin
	# and partial count, fits; otherwise an empty array sends them to the
	# int64 path.
	if ys <= 2147483647:
		counts32 = numpy.empty(nt, dtype=numpy.int32)
		for j in range(nt):
			counts32[j] = rr_counts[j]
	else:
		counts32 = numpy.empty(0, dtype=numpy.int32)
	w = numpy.empty(nt, dtype=numpy.float64)
	for j in range(nt):
		w[j] = rr_counts[j] / ys

	# The binned stage (`gamma_int` column and `f` row) of the most frequent
	# cached columns, per (i_min, bin_scale) class, filled as queries reach
	# them. One copy per thread, so nothing is shared under prange; at most
	# 8 classes and the 4 most frequent columns; the sweep in iteration 65
	# found 4, 16 and 53 columns equally fast.
	n_h = 8
	h_row = 2 * nt + 8 * (n_score_bins + 1)
	n_hs = min(len(q_cached), 4, 2 ** 24 // (n_threads * n_h * h_row))
	_H_keys = numpy.zeros((n_threads, n_h+1, 2), dtype='int64')
	_H_filled = numpy.zeros((n_threads, n_h, n_hs), dtype=numpy.bool_)
	_H_int = numpy.empty((n_threads, n_h, n_hs, nt), dtype='int16')
	_H_f = numpy.empty((n_threads, n_h, n_hs, n_score_bins+1), dtype='float64')

	# Queries are sorted by width (`order`, a stable argsort of Q_lens made
	# by `tomtom` with numpy, which saves compiling a sort), and consecutive
	# queries of one width form batches: the largest power of two, at most
	# n_batch, of those left in the run. Each query's arithmetic and output
	# row are unchanged.
	n_q = len(Q_lens)
	b_start = numpy.empty(n_q+1, dtype='int64')
	n_b = 0
	u = 0
	while u < n_q:
		nq = Q_lens[order[u]]
		n_run = 0
		while u + n_run < n_q and n_run < n_batch and \
				Q_lens[order[u + n_run]] == nq:
			n_run += 1
		n_take = 1
		while n_take * 2 <= n_run:
			n_take *= 2
		b_start[n_b] = u
		n_b += 1
		u += n_take
	b_start[n_b] = n_q

	# Batches are formed before they are given to threads, so a width's
	# queries fill 16-query batches at any thread count; dealing queries to
	# threads first split each run n_threads ways, and at 8 threads 94
	# queries ran alone. With several threads, `tomtom` sets the parallel
	# chunksize to 1 so the scheduler balances batches dynamically; the
	# default, one contiguous block per thread, was 20% slower at 8 threads.
	for bi in prange(n_b):
		pid = numba.get_thread_id()
		qb = numpy.empty(n_batch, dtype='int64')
		offs = numpy.empty(n_batch, dtype='uint64')
		act = numpy.empty(n_batch, dtype=numpy.bool_)
		b_lo = numpy.empty(n_batch, dtype='int64')
		b_col = numpy.empty(n_batch, dtype='int64')

		u0 = b_start[bi]
		nb = b_start[bi+1] - u0
		for b in range(nb):
			qb[b] = order[u0 + b]
		nq = Q_lens[qb[0]]

		G3 = _G[pid, :nt*nq*nb].reshape((nt, nq, nb))

		for b in range(nb):
			i = qb[b]
			g1 = _gamma_int[pid, :nt*nq].reshape((nt, nq))
			# The distances and medians do not depend on the gamma_int
			# type; the binned stage goes into the batch's buffer (`_Gs` or
			# G3) when the offset fits it, which includes a query running
			# alone (it is then copied into g1), and into this query's int16
			# g1 otherwise. The offset is in [0, n_score_bins] (see `tomtom`),
			# so the int16 call is a guard, never taken on any input tried:
			# it runs in object mode and compiles on first use
			# (`_integer_histogram_int16`), and only the batch buffer's
			# specialization is compiled with `_tomtom`.
			i_min, bin_scale, off_s = _distances_and_medians(Q, T, 
				_gamma[pid], _medians[pid], _median_bins[pid], Q_norm, T_norm,
				rr_counts, Q_offsets[i], nq, n_score_bins, q_slot, S_cache,
				halfway, counts32)
			in_g3 = g_max < 0 or _fits_gamma(off_s, n_score_bins, g_max)
			alone = nb == 1 or not in_g3
			if not in_g3:
				_integer_histogram_int16(T, _gamma[pid], g1, _f[pid], 
					_medians[pid], rr_counts, Q_offsets[i], nq, i_min, 
					bin_scale, off_s, q_slot, G_cache, _H_keys[pid], 
					_H_filled[pid], _H_int[pid], _H_f[pid], w)
			else:
				# With an int8 _G, `_Gs` holds each query's packed gamma_int
				# until `_interleave` copies the batch into _G, and `direct`
				# is None. Otherwise `_Gs` is None and the query writes its
				# lane of _G directly. numba prunes a branch on `x is not
				# None` only when x is None, so each call compiles one of the
				# two.
				if _Gs is not None:
					_integer_histogram(T, _gamma[pid], 
						_Gs[pid, b, :nt*nq].reshape((nt, nq)), _f[pid], 
						_medians[pid], rr_counts, Q_offsets[i], nq, i_min, 
						bin_scale, off_s, q_slot, G_cache, _H_keys[pid], 
						_H_filled[pid], _H_int[pid], _H_f[pid], w)
				if direct is not None:
					_integer_histogram(T, _gamma[pid], G3[:, :, b], _f[pid], 
						_medians[pid], rr_counts, Q_offsets[i], nq, i_min, 
						bin_scale, off_s, q_slot, G_cache, _H_keys[pid], 
						_H_filled[pid], _H_int[pid], _H_f[pid], w)
			offset = uint64(off_s)
			offs[b] = offset

			# The backgrounds span nq*(n_score_bins+offset) bins. When the
			# offset exceeds `n_cache` this can overrun the shared
			# workspace, so allocate a large enough one for this query.
			n_needed = nq*n_score_bins + nq*offset
			# A batched query's sums must fit the batch's dtype: int32 as
			# before for an int16 gamma_int, int16 for an int8 one.
			if g_max >= 0:
				fits = _sums_fit(nq, offset, g_max, 32767)
			else:
				fits = int64(nq) * (int64(offset) + 32768) <= 2147483647
			if n_needed > n_len:
				A = numpy.empty((nq, nq, n_needed), dtype='float64')
				B = numpy.empty((T_max+1) * n_needed, dtype='float64')
				A_csum = numpy.empty((nq, nq, n_needed), dtype='float64')
			else:
				n_a = nq*nq*n_needed
				A = _A[pid, :n_a].reshape((nq, nq, n_needed))
				A_csum = _A_csum[pid, :n_a].reshape((nq, nq, n_needed))
				B = _B[pid, b]

			# A_csum's memory holds the packed span rows; A[i, i] is not set.
			# B's rows are stored flat in the layout (b_lo, b_col).
			b_lo[b], b_col[b] = _p_value_backgrounds_windowed(_f[pid], A, B, 
				A_csum, nq, n_score_bins, T_max, offset, needed, 
				A_csum.reshape(A_csum.size))

			# A single query, or one the batch cannot hold, runs alone
			# from a packed copy of its gamma_int (already in g1 when it did
			# not fit the batch buffer).
			act[b] = not alone and n_needed <= n_len and fits
			if not act[b]:
				if in_g3:
					if _Gs is not None:
						Gb = _Gs[pid, b, :nt*nq].reshape((nt, nq))
						for ri in range(nt):
							for li in range(nq):
								g1[ri, li] = Gb[ri, li]
					if direct is not None:
						for ri in range(nt):
							for li in range(nq):
								g1[ri, li] = G3[ri, li, b]
				_p_values(g1, B[:(T_max+1)*b_col[b]].reshape((T_max+1, 
					b_col[b])), rr_inv, T_lens, -1, nq, offset, 
					_results[pid, b], reverse_complement, b_lo[b])

		if nb > 1:
			if _Gs is not None:
				_interleave(_Gs[pid], _G[pid, :nt*nq*nb], nt*nq, nb)
			G2 = _G[pid, :nt*nq*nb].reshape((nt, nq*nb))
			_p_values_batch(G2, _B[pid], rr_inv, T_lens, nq, offs, act,
				_results[pid], nb, s_proto, b_lo, b_col)

		for b in range(nb):
			i = qb[b]
			res = _results[pid, b]

			# The full-matrix, two-strand case merges straight into the
			# output.
			if reverse_complement == 1 and n_nearest == -1:
				_merge_rc_results_into(res, results[i])
			else:
				if reverse_complement == 1:
					_merge_rc_results(res)
				else:
					res[:, 4] = 0

				if n_nearest == -1:
					for ti in range(n_in_targets):
						for ci in range(5):
							results[i, ti, ci] = res[ti, ci]
				else:
					idxs = numpy.argsort(res[:n_in_targets, 0])[:n_nearest]
					for ti in range(len(idxs)):
						for ci in range(5):
							results[i, ti, ci] = res[idxs[ti], ci]
						results[i, ti, 5] = idxs[ti]


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

	# `_tomtom` is compiled once per layout of Q, and Q's layout follows the
	# inputs': (4, w) motifs from `read_meme` give an F-ordered Q, C-ordered
	# motifs a C-ordered one, and each layout costs a full compile (~30 s
	# cold). Q is passed F-ordered always; Q_norm is computed above from the
	# array as given, so no value changes. A no-op for F-ordered inputs.
	Q = numpy.asfortranarray(Q)
	
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

		# Advanced indexing along the last axis returns an F-ordered array;
		# a C-ordered T makes each letter's row contiguous for the blocked
		# distance loops. The values are unchanged.
		T = numpy.ascontiguousarray(T[:, rr_idxs])
		T_norm = T_norm[rr_idxs]
		rr_inv = rr_inv.astype('uint64')
	else:
		# The same layout and dtypes as the hashed branch, so both reuse one
		# compiled `_tomtom` (an F-ordered T and an int64 `rr_inv` compiled
		# a second one, ~33 s cold). T_norm is computed above.
		T = numpy.ascontiguousarray(T)
		rr_inv = numpy.arange(T.shape[-1], dtype='uint64')
		rr_counts = numpy.ones(T.shape[-1], dtype='int64')
	
	# Query columns with identical bytes give identical distance rows, minima,
	# maxima and medians, so a column that occurs more than once has those
	# computed once (`_fill_column_cache`). The cache holds one float64 row of
	# T.shape[-1] per column, so it is capped at 32 MB, filled with the most
	# frequent columns first. `q_cached` lists the first occurrence of each
	# cached column and q_slot[c] is column c's row in the cache, or -1.
	Qc = numpy.ascontiguousarray(Q.T)
	Qv = Qc.view(numpy.dtype((numpy.void, Qc.dtype.itemsize * Qc.shape[1])))
	_, q_first, q_inv, q_counts = numpy.unique(Qv.ravel(), return_index=True,
		return_inverse=True, return_counts=True)
	q_inv = q_inv.ravel()
	n_slots = max(1, 2 ** 25 // (8 * T.shape[-1]))
	shared = numpy.argsort(-q_counts, kind='stable')[:n_slots]
	shared = shared[q_counts[shared] > 1]
	unique_slot = numpy.full(len(q_counts), -1, dtype='int64')
	unique_slot[shared] = numpy.arange(len(shared))
	q_slot = unique_slot[q_inv]
	q_cached = q_first[shared].astype('int64')

	# With n_score_bins <= 127 the batched `gamma_int` is int8 and its sums
	# int16: each value is x - offset with x in [0, n_score_bins], and the
	# offset is in [0, n_score_bins] (`_fits_gamma` checks it per query).
	# Otherwise int16 and int32, as before.
	narrow = n_score_bins <= 127
	s_proto = numpy.empty(0, dtype='int16' if narrow else 'int32')

	# The large per-call arrays are allocated here rather than in `_tomtom`.
	# numpy asks the kernel for transparent huge pages on allocations of
	# 4 MB or more (madvise), numba's allocator does not, so the first touch
	# of these ~400 MB costs 2 MB page faults instead of 4 kB ones: 85k minor
	# faults per call fell to 5k on JASPAR (iteration 95). Shapes are those
	# `_tomtom` used; `n_batch` there is 16.
	n_in = len(T_lens) // 2 if reverse_complement else len(T_lens)
	Q_max, T_len_max, nt = int(Q_lens.max()), int(T_lens.max()), T.shape[-1]
	n_len = Q_max*n_score_bins + Q_max*n_cache
	results = numpy.empty((len(Q_lens), n_in if n_nearest == -1 else
		n_nearest, 5 if n_nearest == -1 else 6), dtype='float64')
	_A = numpy.empty((n_jobs, Q_max*Q_max*n_len), dtype='float64')
	_A_csum = numpy.empty((n_jobs, Q_max*Q_max*n_len), dtype='float64')
	_B = numpy.empty((n_jobs, 16, (T_len_max+1)*n_len), dtype='float64')
	_G = numpy.empty((n_jobs, nt*Q_max*16), dtype='int8' if narrow else 'int16')
	G_cache = numpy.empty((len(q_cached), nt), dtype='float64')
	# With an int8 _G, each batched query's packed `gamma_int` before
	# `_interleave`; rows padded to a multiple of 8 bytes so each starts on
	# an 8-byte boundary. For int16 staging measured slower (iteration 109).
	_Gs = numpy.empty((n_jobs, 16, -(-nt*Q_max // 8) * 8), dtype='int8') \
		if narrow else None

	# With several threads, `_tomtom`'s parallel loop hands out one batch
	# at a time rather than one contiguous block per thread. Set here, not
	# inside `_tomtom`, because calling it there stops numba caching it.
	_chunk = numba.set_parallel_chunksize(1 if n_jobs > 1 else 0)
	try:
		results = _tomtom(Q, T, Q_lens, T_lens, Q_norm, T_norm, rr_inv,
			rr_counts, n_nearest, n_score_bins, n_median_bins, n_cache,
			n_jobs, int(reverse_complement), q_slot, q_cached, results,
			_A, _A_csum, _B, _G, _Gs, None if narrow else True, G_cache,
			127 if narrow else -1, s_proto, numpy.argsort(Q_lens, kind='stable'))
	finally:
		numba.set_parallel_chunksize(_chunk)

	if n_jobs != -1:
		numba.set_num_threads(_n_jobs)

	return results.transpose(2, 0, 1)
