# fimo.py
# Author: Jacob Schreiber <jmschreiber91@gmail.com>

import math
import numba
import numpy
import pandas
import pyfaidx
import time

from .io import read_meme

from tqdm import tqdm

@numba.njit('float64(float64, float64)', cache=True)
def logaddexp2(x, y):
	"""Calculate the logaddexp in a numerically stable manner in base 2.

	This function is a fast implementation of the logaddexp2 function that
	operates on two numbers and is numerically stable. It should mimic the
	functionality of numpy.logaddexp2 except that it does not have the overhead
	of working on numpy arrays.


	Parameters
	----------
	x: float32
		A single number in log space.

	y: float32
		Another single number in log space.


	Returns
	-------
	z: float32
		The result of log2(pow(2, x) + pow(2, y))
	"""

	if x == float("-inf") and y == float("-inf"):
		return float("-inf")

	if x == float("inf") or y == float("inf"):
		return float("inf")

	vmax, vmin = max(x, y), min(x, y)
	return vmax + math.log2(math.pow(2, vmin - vmax) + 1)


@numba.njit(cache=True)
def _pwm_to_mapping(log_pwm, bin_size):
	"""An internal method for calculating score <-> log p-value mappings.

	This function takes in a PWM consisting of log probabilities and outputs
	a mapping between observed scores (as a convolution of the PWM across a
	one-hot encoded sequence) and log p-values. This mapping is calculated 
	quickly using dynamic programming scanning over all potential sequences.

	Importantly, the p-values are in log space meaning that values near zero
	at the start of the array are insignificant whereas those with large 
	magnitude towards the end of the array are more statistically significant.


	Parameters
	----------
	log_pwm: numpy.ndarray, shape=(len(alphabet), length)
		A position-weight matrix containing a motif encoded as the log
		probability of any character in any position.

	bin_size: float
		The size of the score bins to map to p-values. The smaller this value,
		the more bins, indicating higher precision but also longer calculation
		time.


	Returns
	-------
	smallest: int
		The number of bins between true zero and the smallest value in the
		array. In other words, the offset to subtract from binned scores to get
		p-values.

	log1mcdf: numpy.ndarray
		The log of 1 minus the cdf, or in other words, the log p-values
		associated with each score bin.
	"""

	int_log_pwm = numpy.round(log_pwm / bin_size).astype(numpy.int32)
	smallest, largest = _mapping_range(int_log_pwm)

	logpdf = numpy.empty(largest - smallest + 1)
	old_logpdf = numpy.empty(largest - smallest + 1)
	_fill_mapping(int_log_pwm, smallest, logpdf, old_logpdf)
	return smallest, old_logpdf


@numba.njit(cache=True)
def _mapping_range(int_log_pwm):
	"""The smallest and largest score bins of `_pwm_to_mapping`'s table for a
	PWM already divided by the bin size, rounded and cast to int32. The table
	has largest - smallest + 1 bins."""

	n, l = int_log_pwm.shape

	smallest, largest = 9999999, -9999999
	log_pwm_min_csum, log_pwm_max_csum = 0, 0
	for i in range(l):
		log_pwm_min = 9999999
		log_pwm_max = -9999999

		for j in range(n):
			log_pwm_min = min(log_pwm_min, int_log_pwm[j, i])
			log_pwm_max = max(log_pwm_max, int_log_pwm[j, i])

		log_pwm_min_csum += log_pwm_min
		log_pwm_max_csum += log_pwm_max

		smallest = min(smallest, log_pwm_min_csum)
		largest = max(largest, log_pwm_max_csum)

	largest += l
	return smallest, largest


@numba.njit(cache=True)
def _fill_mapping(int_log_pwm, smallest, logpdf, old_logpdf):
	"""The dynamic program of `_pwm_to_mapping`, leaving the table of log
	p-values in `old_logpdf`. `logpdf` is scratch. Both have one entry per bin
	from `_mapping_range` and may hold anything on entry.

	Every value in the table is -inf or finite, and never -0.0: a sum is -0.0
	only when both terms are, and each value is either log_bg or a logaddexp2
	of two such values, `vmax + log2(2 ** (vmin - vmax) + 1)` with a second
	term of at least +0.0. So `logaddexp2(-inf, y)` is `y + log2(pow(2, -inf) +
	1) = y + 0.0 = y` bitwise, as `pow(2, -inf)` is +0.0 and `log2(1)` is +0.0
	(C99 Annex F), and in the same way `logaddexp2(x, -inf)` is `x`. Those
	calls are written as their result; every other call is made as before, on
	the same bins in the same order.
	"""

	n, l = int_log_pwm.shape
	log_bg = math.log2(0.25)
	size = old_logpdf.shape[0]

	for j in range(size):
		old_logpdf[j] = -numpy.inf

	for i in range(n):
		idx = int_log_pwm[i, 0] - smallest
		if old_logpdf[idx] == -numpy.inf:
			old_logpdf[idx] = log_bg
		else:
			old_logpdf[idx] = logaddexp2(old_logpdf[idx], log_bg)

	for i in range(1, l):
		for j in range(size):
			logpdf[j] = -numpy.inf

		for j, x in enumerate(old_logpdf):
			if x != -numpy.inf:
				y = log_bg + x
				for k in range(n):
					idx = j + int_log_pwm[k, i]
					if logpdf[idx] == -numpy.inf:
						logpdf[idx] = y
					else:
						logpdf[idx] = logaddexp2(logpdf[idx], y)

		for j in range(size):
			old_logpdf[j] = logpdf[j]

	# The survival function, summed from the top bin down.
	for i in range(size - 2, -1, -1):
		if old_logpdf[i] == -numpy.inf:
			old_logpdf[i] = old_logpdf[i + 1]
		elif old_logpdf[i + 1] != -numpy.inf:
			old_logpdf[i] = logaddexp2(old_logpdf[i], old_logpdf[i + 1])


@numba.njit(cache=True)
def _table_layout(motifs, motif_lengths, bin_size):
	"""Where each motif's p-value table goes in one flat array, for
	`_fill_tables`.

	Returns the PWMs divided by the bin size, rounded and cast to int32 as
	`_pwm_to_mapping` does (the same operations on each element), each motif's
	smallest bin, the offsets of its table, and the cost of computing it:
	about its width times its number of bins. A motif with no columns gets an
	empty table.
	"""

	n = len(motif_lengths) - 1
	int_pwms = numpy.round(motifs / bin_size).astype(numpy.int32)

	smallests = numpy.zeros(n, dtype=numpy.int64)
	offsets = numpy.empty(n + 1, dtype=numpy.int64)
	cost = numpy.empty(n, dtype=numpy.int64)
	offsets[0] = 0
	for i in range(n):
		s, e = motif_lengths[i], motif_lengths[i+1]
		size = 0
		if e > s:
			smallest, largest = _mapping_range(int_pwms[:, s:e])
			smallests[i] = smallest
			size = largest - smallest + 1

		offsets[i+1] = offsets[i] + size
		cost[i] = (numpy.int64(e) - numpy.int64(s)) * size

	return int_pwms, smallests, offsets, cost


@numba.njit(parallel=True, cache=True)
def _fill_tables(int_pwms, motif_lengths, smallests, offsets, order, tables):
	"""Write the table of each motif, as `_pwm_to_mapping` returns it, into
	`tables[offsets[i]:offsets[i+1]]`."""

	n = len(motif_lengths) - 1
	for t in numba.prange(n):
		i = order[t]
		s, e = motif_lengths[i], motif_lengths[i+1]
		a, b = offsets[i], offsets[i+1]
		if b > a:
			logpdf = numpy.empty(b - a)
			_fill_mapping(int_pwms[:, s:e], smallests[i], logpdf, tables[a:b])


def _pvalue_tables(motifs, motif_lengths, bin_size):
	"""The p-value tables of `_all_pwm_to_mapping`, concatenated: each motif's
	smallest bin, the offsets of its table, and the tables in one array."""

	int_pwms, smallests, offsets, cost = _table_layout(motifs, motif_lengths,
		bin_size)

	# A prange hands each thread one contiguous block of iterations, so the
	# motifs are sorted by cost and dealt round-robin into one block per
	# thread: block t is by_cost[t], by_cost[t + n_blocks], ... The order
	# changes no value. It is built here because numba's argsort adds about
	# 0.6 s to a cold compile.
	n = len(cost)
	n_blocks = max(min(numba.get_num_threads(), n), 1)
	by_cost = numpy.full(-(-n // n_blocks) * n_blocks, -1, dtype=numpy.int64)
	by_cost[:n] = numpy.argsort(-cost, kind='stable')
	order = by_cost.reshape(-1, n_blocks).T.flatten()
	order = order[order >= 0]
	tables = numpy.empty(offsets[-1], dtype=numpy.float64)
	_fill_tables(int_pwms, motif_lengths, smallests, offsets, order, tables)
	return smallests, offsets, tables


@numba.njit(parallel=True, cache=True)
def _all_pwm_to_mapping(motifs, motif_lengths, bin_size):
	n = len(motif_lengths) - 1

	smallests = numpy.empty(n, dtype='int64')
	logpdfs = [numpy.empty(0) for i in range(n)]

	for i in numba.prange(n):
		s, e = motif_lengths[i], motif_lengths[i+1]

		smallest, logpdf = _pwm_to_mapping(motifs[:, s:e], bin_size)
		smallests[i] = smallest
		logpdfs[i] = logpdf

	return smallests, logpdfs


# How `_fast_hits` bounds a window. Every sequence position gets a code for
# the q letters starting there, and each motif gets tables of the summed
# weights of up to q consecutive columns for every code, so a block of
# columns costs one lookup. A code exists at every position, so a block can
# start at any column. q is the largest value up to `_QMAX` whose tables have
# at most `_TABLE_MAX` entries: 5 for DNA, whose 4 letters plus N give
# 5**5 = 3125. Codes are uint16, so `_TABLE_MAX` must stay at most 65536.
#
# `_block_layout` splits each motif's columns into blocks. A window is first
# tested on two blocks, for a motif of at least 2q columns the two full blocks
# with the largest summed gap between each column's largest and mean entry,
# and then after every further block, in the same order of gap. Chosen by
# measurement on the benchmark: the first test on the best two blocks passes
# 0.27% of windows against 1.02% for the motif's first 2q columns. A first
# test of one lookup passes three times as many as the first 2q columns, and
# its loop was no faster per window than the loop of two.
_QMAX = 5
_TABLE_MAX = 3125


def _score_bounds(pwm, pwm_lengths, thresholds):
	"""Upper bounds that let `_fast_hits` abandon a window early.

	A column can add at most its largest entry to a window's score, or 0.0 when
	the position is an N, which reads the kernel's all-zero row: that is
	`col_max`, per global column, with +inf for a column holding NaN.
	`_block_layout` sums it over the blocks a window has not yet been tested
	on. `gap` is each column's largest entry minus its mean, which only
	decides the order of the blocks. `tops[k]` is the largest score motif k can
	reach at all.

	A window is abandoned only when its partial score plus the remaining bound
	is at most `cuts[k]`, the threshold minus a margin. The partial score is a
	sum of table entries that each sum up to q columns (with +0.0 for columns
	outside the block, which is exact), and the blocks tested and the blocks
	left partition the motif's columns, whatever order they are visited in.
	So the partial score plus the remaining bound is one more order of addition
	of one term per column, each at least that column's entry. Every float sum
	involved (the partial score, the bounds, and the full score, in any order
	of addition) has at most n terms of magnitude at most
	W = sum over columns of the largest finite |entry|, and any order of
	summation is within (n - 1) * 2**-53 * W of the exact sum (to first order).
	Chaining these through the test, plus the rounding of the test's own
	addition and of the cut, puts the full score of an abandoned window at most
	(3n + 1) * 2**-53 * (W + |threshold|) above threshold - margin. The margin,
	(n + 2) * 1e-12 * (1 + W + |threshold|), is over 3,000 times that, and far
	below the smallest gap between a window's score and its threshold on real
	motifs (4.7e-5 on the benchmark).

	An unreachable threshold (+inf) gets a cut of +inf, so every window is
	abandoned, as `_fast_hits` would report none of them. Any other non-finite
	cut, and a motif with no columns, gets -inf: nothing is abandoned and the
	full score decides, as before.
	"""

	lengths = numpy.asarray(pwm_lengths, dtype=numpy.int64)
	widths = numpy.diff(lengths)
	n_motifs = len(widths)
	max_width = int(widths.max()) if n_motifs > 0 else 0
	mask = numpy.arange(max_width) < widths[:, None]

	with numpy.errstate(invalid='ignore'):
		raw_max = pwm.max(axis=0)
		col_max = numpy.where(numpy.isnan(raw_max), numpy.inf,
			numpy.maximum(raw_max, 0.0)).astype(numpy.float64)
		col_abs = numpy.where(numpy.isfinite(pwm), numpy.abs(pwm), 0.0).max(
			axis=0)
		gap = (raw_max - pwm.mean(axis=0)).astype(numpy.float64)

	# One row per motif, zero-padded, so each sum is over its own motif.
	padded = numpy.zeros((n_motifs, max_width + 1))
	padded[:, :max_width][mask] = col_max
	suffix = numpy.cumsum(padded[:, ::-1], axis=1)[:, ::-1]
	tops = numpy.where(widths > 0, suffix[:, 0], numpy.inf)

	padded = numpy.zeros((n_motifs, max_width))
	padded[mask] = col_abs
	magnitude = padded.sum(axis=1)

	t = numpy.asarray(thresholds, dtype=numpy.float64)
	with numpy.errstate(invalid='ignore', over='ignore'):
		margin = (widths + 2) * 1e-12 * (1.0 + magnitude + numpy.abs(t))
		cuts = t - margin
	cuts = numpy.where(numpy.isfinite(cuts) & (widths > 0), cuts, -numpy.inf)
	cuts = numpy.where((t == numpy.inf) & (widths > 0), numpy.inf, cuts)
	return col_max, gap, cuts, tops


def _qmer_width(n_rows):
	"""The number of letters per code, for a PWM with `n_rows` rows including
	the zero row for N: the largest q up to `_QMAX` with n_rows**q entries at
	most `_TABLE_MAX`."""

	q = 1
	while q < _QMAX and n_rows ** (q + 1) <= _TABLE_MAX:
		q += 1
	return q


@numba.njit(cache=True)
def _qmer_codes(X, n_rows, q):
	"""The code of the q letters starting at each position of `X`,
	X[i] + n_rows * X[i+1] + ... + n_rows**(q-1) * X[i+q-1]. Positions past
	the end of `X` read as N (n_rows - 1), and there are 2q codes more than
	positions, so a read at a window's start plus q is always in bounds. A code
	that runs past the end of its sequence is only ever read through table
	columns that weigh 0.0. q = 5 and 4 are written out, which vectorizes;
	another q takes the same arithmetic in the loop at the end."""

	L = X.shape[0]
	n_codes = L + 2 * q
	codes = numpy.empty(n_codes, dtype=numpy.uint16)
	b = numpy.uint16(n_rows)
	letter_n = numpy.uint16(n_rows - 1)

	m = max(L - q + 1, 0)
	if q == 5:
		for i in range(m):
			codes[i] = (((numpy.uint16(X[i+4]) * b + numpy.uint16(X[i+3])) * b
				+ numpy.uint16(X[i+2])) * b + numpy.uint16(X[i+1])) * b + \
				numpy.uint16(X[i])
	elif q == 4:
		for i in range(m):
			codes[i] = ((numpy.uint16(X[i+3]) * b + numpy.uint16(X[i+2])) * b
				+ numpy.uint16(X[i+1])) * b + numpy.uint16(X[i])
	else:
		m = 0

	for i in range(m, n_codes):
		c = numpy.uint16(0)
		for d in range(q - 1, -1, -1):
			if i + d < L:
				c = c * b + numpy.uint16(X[i + d])
			else:
				c = c * b + letter_n
		codes[i] = c

	return codes


@numba.njit(cache=True)
def _block_layout(col_max, gap, pwm_lengths, q):
	"""Split each motif's columns into blocks of at most q consecutive columns,
	in the order `_fast_hits` tests them.

	Block h of motif k covers columns `blk_off[h] .. blk_off[h] + blk_w[h] - 1`
	of the motif, for h from `blk_start[k]` to `blk_start[k+1] - 1`, and the
	blocks partition its columns. The first test reads the first two blocks.
	For a motif of at most q columns they are the whole motif and an empty
	block at column 0, whose table is all 0.0; for one of fewer than 2q columns
	columns 0 to q-1 and q onwards; and for a wider motif the pair of
	non-overlapping full blocks with the largest summed `gap`, the first of
	equals. The columns left over are cut into blocks of q from the left of
	each run of consecutive columns, and tested largest summed gap first,
	keeping their order among equals. A NaN gap compares as never larger, so
	it only changes which valid layout is used.

	`blk_rest[h]` is the sum of `col_max` over the blocks after h: the most
	those columns can still add to a window's score. It is 0.0 for each motif's
	last block. Any order of summation is covered by `_score_bounds`' margin.
	"""

	n_motifs = len(pwm_lengths) - 1
	q = numpy.int64(q)
	cap = numpy.int64(pwm_lengths[n_motifs]) + 2 * n_motifs
	blk_off = numpy.empty(cap, dtype=numpy.uint64)
	blk_w = numpy.empty(cap, dtype=numpy.uint64)
	blk_rest = numpy.empty(cap, dtype=numpy.float64)
	blk_start = numpy.empty(n_motifs + 1, dtype=numpy.int64)
	key = numpy.empty(cap, dtype=numpy.float64)
	seg_lo = numpy.empty(3, dtype=numpy.int64)
	seg_hi = numpy.empty(3, dtype=numpy.int64)

	b = numpy.int64(0)
	for k in range(n_motifs):
		off = numpy.int64(pwm_lengths[k])
		n = numpy.int64(pwm_lengths[k + 1]) - off
		blk_start[k] = b

		n_seg = 0
		if n <= q:
			blk_off[b], blk_w[b] = 0, n
			blk_off[b + 1], blk_w[b + 1] = 0, 0
		elif n < 2 * q:
			blk_off[b], blk_w[b] = 0, q
			blk_off[b + 1], blk_w[b + 1] = q, n - q
		else:
			# The summed gap of the full block starting at each column.
			n_full = n - q + 1
			gs = numpy.empty(n_full, dtype=numpy.float64)
			for o in range(n_full):
				v = 0.0
				for d in range(q):
					v += gap[off + o + d]
				gs[o] = v

			o1, o2 = 0, q
			best = gs[0] + gs[q]
			for o in range(n_full):
				for p in range(o + q, n_full):
					v = gs[o] + gs[p]
					if v > best:
						best = v
						o1, o2 = o, p

			blk_off[b], blk_w[b] = o1, q
			blk_off[b + 1], blk_w[b + 1] = o2, q
			seg_lo[0], seg_hi[0] = 0, o1
			seg_lo[1], seg_hi[1] = o1 + q, o2
			seg_lo[2], seg_hi[2] = o2 + q, n
			n_seg = 3
		b += 2

		# The columns left over, in blocks of at most q, largest gap first.
		r0 = b
		for e in range(n_seg):
			j = seg_lo[e]
			while j < seg_hi[e]:
				w = min(q, seg_hi[e] - j)
				v = 0.0
				for d in range(w):
					v += gap[off + j + d]
				h = b
				while h > r0 and v > key[h - 1]:
					blk_off[h], blk_w[h], key[h] = blk_off[h-1], blk_w[h-1], key[h-1]
					h -= 1
				blk_off[h], blk_w[h], key[h] = j, w, v
				b += 1
				j += w

		s = 0.0
		for h in range(b - 1, blk_start[k] - 1, -1):
			blk_rest[h] = s
			for d in range(numpy.int64(blk_w[h])):
				s += col_max[off + numpy.int64(blk_off[h]) + d]

	blk_start[n_motifs] = b
	return blk_off[:b].copy(), blk_w[:b].copy(), blk_rest[:b].copy(), blk_start


@numba.njit(inline='always')
def _qmer_tables(tab, pwm, off, q, stride, blk_off, blk_w, b0, n_blocks):
	"""tab[g * stride + c] = the summed weights of the columns of block b0 + g
	of the motif whose columns start at `off`, for the q letters with code c:
	digit d of c is the letter at the block's column d. `stride` is at least
	the number of codes, n_rows**q. Digits at or past the block's width weigh
	0.0. Each table is built in place, from its last digit to its first, each
	digit multiplying its size by the number of rows."""

	n_rows = numpy.int64(pwm.shape[0])
	q = numpy.int64(q)
	off = numpy.int64(off)
	stride = numpy.int64(stride)
	b0 = numpy.int64(b0)
	for g in range(numpy.int64(n_blocks)):
		t0 = g * stride
		col0 = off + numpy.int64(blk_off[b0 + g])
		w = numpy.int64(blk_w[b0 + g])
		tab[t0] = 0.0
		size = numpy.int64(1)
		for d in range(q - 1, -1, -1):
			for t in range(size - 1, -1, -1):
				v = tab[t0 + t]
				for s in range(n_rows - 1, -1, -1):
					wt = 0.0
					if d < w:
						wt = numpy.float64(pwm[s, col0 + d])
					tab[t0 + t * n_rows + s] = v + wt
			size *= n_rows


@numba.njit(inline='always')
def _copy_hits(hits, o, seqs, starts, ends, scores, pvals):
	"""Copy one motif's hits into rows `o` onwards of the hit columns."""

	for h in range(len(hits)):
		l, start, end, score, pval = hits[h]
		seqs[o+h] = l
		starts[o+h] = numpy.int64(start)
		ends[o+h] = numpy.int64(end)
		scores[o+h] = score
		pvals[o+h] = pval


@numba.njit(parallel=True, fastmath=True, cache=True)
def _fast_hits(X, codes, q, chrom_lengths, pwm, pwm_lengths, score_threshold,
	bin_size, smallest, score_to_pvals, score_to_pval_lengths, blk_off, blk_w,
	blk_rest, blk_start, cuts, tops, order):
	"""Scan every motif over every sequence and return the hits as columns.

	`codes` and `q` come from `_qmer_codes`: the bound that abandons a window
	reads a block of up to q columns per table lookup, and the blocks and the
	order they are read in come from `_block_layout`. The hits of each motif
	are collected in a list while it is scanned, then copied into flat arrays
	in the order given by `order`, so that the hits of motif `order[t]` are
	rows `offsets[t]` to `offsets[t+1]` of the returned sequence index, start,
	end, score and p-value arrays.
	"""

	n_motifs = len(pwm_lengths) - 1
	n_chroms = len(chrom_lengths) - 1
	q = numpy.uint64(q)

	hits = []
	for i in range(n_motifs):
		j = numpy.int64(1)
		k = numpy.uint64(1)
		l = numpy.float64(1.0)
		hits.append([(j, k, k, l, l) for z in range(0)])

	for k in numba.prange(n_motifs):
		n = pwm_lengths[k+1] - pwm_lengths[k]
		b0 = numpy.uint64(blk_start[k])
		n_blocks = numpy.uint64(blk_start[k+1] - blk_start[k])
		k = numpy.uint64(k)
		thresh = score_threshold[k]
		cut = cuts[k]
		off = numpy.uint64(pwm_lengths[k])

		# Every block's table starts `_TABLE_MAX` entries after the last,
		# whatever the number of codes, so the second lookup of the first test
		# is a constant offset from the table's start and needs no register of
		# its own. It is set inside the prange body because a value set outside
		# is passed into the body as an argument, which LLVM cannot fold.
		stride = numpy.uint64(_TABLE_MAX)

		# Skip a motif whose best possible score cannot pass its threshold.
		if tops[k] > cut:
			# One table per block. An N reads the all-zero last row of `pwm`
			# and adds +0.0, in the tables and in the rescoring loop below.
			tab = numpy.empty(n_blocks * stride, dtype=numpy.float64)
			_qmer_tables(tab, pwm, off, q, stride, blk_off, blk_w, b0, n_blocks)

			# The first test reads blocks b0 and b0 + 1, at columns o1 and o2
			# of the window. A block lies inside the motif, so its code is
			# read inside the window.
			rest_p = blk_rest[b0 + numpy.uint64(1)]
			o1 = blk_off[b0]
			o2 = blk_off[b0 + numpy.uint64(1)]

			for l in range(n_chroms):
				start = numpy.uint64(chrom_lengths[l])
				end = numpy.uint64(chrom_lengths[l+1])

				# Windows start at positions x from `start` to `xe` - 1.
				m = numpy.uint64(max(numpy.int64(end - start) -
					numpy.int64(n) + 1, 0))
				xe = start + m
				x = start
				while x < xe:
					# The first test, in a loop of its own that moves to the
					# next window until one survives it, which about 0.3% of
					# windows do on real motifs. With nothing else in the loop,
					# its pointers, position and bounds stay in registers; in
					# one loop with the code below, LLVM reloaded and spilled
					# them for every window. Only the bound test uses these
					# columns.
					bound = 0.0
					while x < xe:
						bound = tab[numpy.uint64(codes[x+o1])] + tab[stride +
							numpy.uint64(codes[x+o2])]
						if bound + rest_p <= cut:
							x += numpy.uint64(1)
							continue
						break

					if x >= xe:
						break

					base = x
					i = x - start
					alive = True
					g = numpy.uint64(2)
					while g < n_blocks:
						bound += tab[g * stride + numpy.uint64(codes[base +
							blk_off[b0+g]])]
						if bound + blk_rest[b0+g] <= cut:
							alive = False
							break
						g += numpy.uint64(1)

					if not alive:
						x += numpy.uint64(1)
						continue

					# A window that survives is scored strictly left to right,
					# as it would be without the bound test. `score` starts at
					# +0.0 and a sum is -0.0 only when both terms are, so it is
					# never -0.0 and adding the +0.0 of an N leaves it bitwise
					# unchanged: the same as skipping the column.
					score = 0.0
					for j in range(n):
						j = numpy.uint64(j)
						idx = numpy.uint64(X[start+i+j])
						m_idx = numpy.uint64(j + pwm_lengths[k])
						score += pwm[idx, m_idx]

					if score > thresh:
						score_idx = int(score / bin_size) - smallest[k]                    
						score_idx += score_to_pval_lengths[k]
						hits[k].append((numpy.int64(l), i, i+n, score, 
							2.0 ** score_to_pvals[score_idx]))

					x += numpy.uint64(1)

	# `numpy.zeros` and a prange each start a parallel region, which costs
	# about 20 us, so the offsets are filled serially and a small number of
	# hits is copied serially.
	offsets = numpy.empty(n_motifs + 1, dtype=numpy.int64)
	offsets[0] = 0
	for t in range(n_motifs):
		offsets[t+1] = offsets[t] + len(hits[order[t]])

	n_total = offsets[n_motifs]
	seqs = numpy.empty(n_total, dtype=numpy.int64)
	starts = numpy.empty(n_total, dtype=numpy.int64)
	ends = numpy.empty(n_total, dtype=numpy.int64)
	scores = numpy.empty(n_total, dtype=numpy.float64)
	pvals = numpy.empty(n_total, dtype=numpy.float64)

	if n_total < 1024:
		for t in range(n_motifs):
			_copy_hits(hits[order[t]], offsets[t], seqs, starts, ends, scores, 
				pvals)
	else:
		for t in numba.prange(n_motifs):
			_copy_hits(hits[order[t]], offsets[t], seqs, starts, ends, scores, 
				pvals)

	return offsets, seqs, starts, ends, scores, pvals


@numba.njit(cache=True)
def _fast_convert(X, mapping):
	for i in range(X.shape[0]):
		X[i] = mapping[X[i]]


@numba.njit(cache=True)
def _int_one_hot_to_index(X, out):
	"""`_one_hot_to_index` for signed integers, with four channels.

	numpy computes the expression in int64, wrapping on overflow, and keeps
	the low byte, which depends only on the low byte of each value. Summing
	the values' low bytes gives the same low byte and cannot overflow, which
	matters because numba's signed arithmetic assumes it does not.
	"""

	n, _, l = X.shape
	for i in range(n):
		o = i * l
		for j in range(l):
			v0 = X[i, 0, j]
			v1 = X[i, 1, j]
			v2 = X[i, 2, j]
			v3 = X[i, 3, j]
			s = (numpy.int64(numpy.int8(v0)) + numpy.int8(v1) + numpy.int8(v2)
				+ numpy.int8(v3))

			# numpy's argmax: the first of the largest values.
			best = v0
			arg = 0
			if v1 > best:
				best = v1
				arg = 1
			if v2 > best:
				best = v2
				arg = 2
			if v3 > best:
				arg = 3

			out[o + j] = numpy.int8((arg + 1) * s - 1)


@numba.njit(cache=True)
def _float_one_hot_to_index(X, sums, out):
	"""`_one_hot_to_index` for floats and unsigned integers, with four channels.

	`sums` is numpy's own `X.sum(axis=1)`, so the order of summation is
	numpy's. numpy computes the rest in float64 and casts it to int8, which
	truncates toward zero; a value outside (-129, 128), or NaN, has no exact
	int8 and casts in a platform-dependent way. Returns False when any value
	falls outside, and the caller then uses numpy for the whole array.

	A NaN in any channel makes the sum NaN, so the position falls back and
	numpy's rule for NaN in argmax is never needed here.
	"""

	n, _, l = X.shape
	ok = True
	for i in range(n):
		o = i * l
		for j in range(l):
			v0 = X[i, 0, j]
			v1 = X[i, 1, j]
			v2 = X[i, 2, j]
			v3 = X[i, 3, j]

			# numpy's argmax: the first of the largest values.
			best = v0
			arg = 0
			if v1 > best:
				best = v1
				arg = 1
			if v2 > best:
				best = v2
				arg = 2
			if v3 > best:
				arg = 3

			v = numpy.float64(arg + 1) * numpy.float64(sums[i, j]) - 1.0
			if v > -129.0 and v < 128.0:
				out[o + j] = numpy.int8(numpy.int64(v))
			else:
				ok = False

	return ok


def _one_hot_to_index(sequences):
	"""Convert one-hot sequences to flattened int8 indices and their offsets.

	The indices are exactly
	`(((sequences.argmax(axis=1) + 1) * sequences.sum(axis=1)) - 1)`
	cast to int8 and flattened: the channel with the largest value, or -1
	where no channel is set. Four-channel arrays of a native-endian integer,
	bool, float32 or float64 dtype are converted in one pass with numba;
	anything else, or a float result that has no exact int8, uses the numpy
	expression. On every path the indices are a new array, never a view of
	`sequences`, so `fimo` may modify them in place.
	"""

	dtype = sequences.dtype
	if sequences.ndim == 3 and sequences.shape[1] == 4 and dtype.isnative:
		n, _, l = sequences.shape
		X = numpy.empty(n * l, dtype=numpy.int8)
		X_lengths = (numpy.arange(n + 1) * l).astype(numpy.int64)

		if dtype == numpy.bool_:
			# numpy sums a bool as 0 or 1 in int64, as the int8 view does.
			_int_one_hot_to_index(sequences.view(numpy.int8), X)
			return X, X_lengths
		elif dtype.kind == 'i':
			_int_one_hot_to_index(sequences, X)
			return X, X_lengths
		elif dtype.kind == 'u' or dtype in (numpy.float32, numpy.float64):
			if _float_one_hot_to_index(sequences, sequences.sum(axis=1), X):
				return X, X_lengths

	X = ((sequences.argmax(axis=1) + 1) * sequences.sum(axis=1)) - 1
	X_lengths = numpy.arange(X.shape[0]+1) * X.shape[-1]

	X = X.astype(numpy.int8).flatten()
	X_lengths = X_lengths.astype(numpy.int64)
	return X, X_lengths


def fimo(motifs, sequences, alphabet=['A', 'C', 'G', 'T'], bin_size=0.1, 
	eps=0.0001, threshold=0.0001, reverse_complement=True, return_counts=False, 
	dim=0):
	"""An implementation of the FIMO algorithm from the MEME suite.

	This function implements the "Finding Individual Motif Instances" (FIMO)
	algorithm from the MEME suite. This algorithm takes a set of PWMs and
	identifies where these PWMs have statistically significant hits against a
	set of sequences. These sequences can either come from a FASTA file, such
	as an entire genome or a set of peaks, or can be one-hot encoded sequences.

	This implementation uses numba to accelerate the inner loop, and
	parallelizes across the motif axis. No support exists for calculating
	q-values as, in my opinion, q-values do not make sense here and are both
	compute- and memory-inefficient.


	Parameters
	----------
	motifs: str or dict
		A MEME file to load containing motifs to scan, or a dictionary where
		the keys are names of motifs and the values are PWMs with shape
		(len(alphabet), pwm_length).

	sequences: str or numpy.ndarray
		A set of sequences to scan the motifs against. If this is a string,
		assumes it is a filepath to a FASTA-formatted file. If this is a numpy
		array, will use those instead.

	alphabet: list, optional
		A list of characters to use for the alphabet, defining the order that
		characters should appear. Default is ['A', 'C', 'G', 'T'].

	bin_size: float, optional
		The size of the bins discretizing the PWM scores. The smaller the bin
		size the higher the resolution, but the less data may be available to
		support it. Default is 0.1.

	eps: float, optional
		A small pseudocount to add to the motif PWMs before taking the log.
		Default is 0.0001.

	threshold: float, optional
		The p-value threshold to use for reporting matches. Default is 0.0001.

	reverse_complement: bool, optional
		Whether to scan each motif and also the reverse complements. Default
		is True.

	return_counts: bool, optioal
		Whether to only return the count of the number of matches instead of
		dataframes containing information about each match. If True, the return
		will be a single array. Default is False

	dim: 0 or 1, optional
		Whether to return one dataframe for each motif containing all hits for
		that motif across all examples (0, default) or one dataframe for each 
		example containing all hits across all motifs to that example (1).
		Default is 0.


	Returns
	-------
	hits: list of pandas.DataFrames or numpy.ndarray
		A list of pandas.DataFrames containing motif hits, where the exact
		semantics of each dataframe are determined by `dim`. Alternatively,
		a numpy array of just the number of counts per motif if return_counts
		is set to True.
	"""

	log_threshold = math.log2(threshold)

	# Extract the motifs and potentially the reverse complements
	if isinstance(motifs, str):
		motifs_ = read_meme(motifs)
	elif isinstance(motifs, dict):
		motifs_ = motifs
	else:
		raise ValueError("`motifs` must be a dict or a filename.")

	motifs_fwd: list[tuple[str, numpy.ndarray]] = []
	motifs_rev: list[tuple[str, numpy.ndarray]] = []
	for name in motifs_:
		# If provided, motifs must be dict[str, np.ndarray]
		pwm = motifs_[name]
		if not isinstance(pwm, numpy.ndarray):
			try:
				pwm = pwm.numpy()
			except:
				raise ValueError(
					f"`motifs` must be a dict[str, numpy.ndarray], not {type(pwm)}."
				)
		motifs_fwd.append((name, pwm))
		
		if reverse_complement:
			motifs_rev.append((name + "-rc", pwm[::-1, ::-1]))
	
	if reverse_complement:
		motifs = [*motifs_fwd, *motifs_rev]
	else:
		motifs = motifs_fwd

	# Initialize arrays to store motif properties
	n_motifs = len(motifs)

	motif_names = numpy.array([name for name, _ in motifs])
	motif_lengths = [0] + [pwm.shape[-1] for _, pwm in motifs]
	motif_lengths = numpy.cumsum(motif_lengths).astype(numpy.uint64)

	motif_pwms = numpy.concatenate([pwm for _, pwm in motifs], axis=-1)
	motif_pwms = numpy.log2(motif_pwms + eps) - math.log2(0.25)

	_smallest, _score_to_pvals_lengths, _score_to_pvals = _pvalue_tables(
		motif_pwms, motif_lengths, bin_size)

	# Each motif's score threshold is the first bin of its table whose log
	# p-value is below the threshold, or inf when no bin is. The first
	# passing bin at or after each table's start is found for all motifs at
	# once, and it belongs to the motif only if it lies before the table's end.
	starts = _score_to_pvals_lengths[:-1]
	passing = numpy.flatnonzero(_score_to_pvals < log_threshold)
	pos = numpy.searchsorted(passing, starts)
	if len(passing) > 0:
		first = passing[numpy.minimum(pos, len(passing) - 1)]
	else:
		first = starts

	found = (pos < len(passing)) & (first < _score_to_pvals_lengths[1:])
	_score_thresholds = numpy.full(n_motifs, numpy.inf, dtype=numpy.float32)
	_score_thresholds[found] = ((first - starts + _smallest) * bin_size)[found]

	# Extract the sequence from a FASTA
	if isinstance(sequences, str):
		fasta = pyfaidx.Fasta(sequences)
		sequence_names = numpy.array(list(fasta.keys()))
		X, lengths = [], [0]
		
		alphabet = ''.join(alphabet)
		alpha_idxs = numpy.frombuffer(bytearray(alphabet, 'utf8'), 
			dtype=numpy.int8)
		one_hot_mapping = numpy.zeros(256, dtype=numpy.int8) - 1
		for i, idx in enumerate(alpha_idxs):
			one_hot_mapping[idx] = i
		
		for name, chrom in fasta.items():
			chrom = chrom[:].seq.upper()
			lengths.append(lengths[-1] + len(chrom))
			
			X_idxs = numpy.frombuffer(bytearray(chrom, "utf8"), 
				dtype=numpy.int8)
			_fast_convert(X_idxs, one_hot_mapping)
			X.append(X_idxs)
			
		X = numpy.concatenate(X)
		X_lengths = numpy.array(lengths, dtype=numpy.int64)
		
	elif not isinstance(sequences, numpy.ndarray):
		sequences = sequences.numpy()
			
	if isinstance(sequences, numpy.ndarray):
		sequence_names = None
		X, X_lengths = _one_hot_to_index(sequences)

	# The kernel's PWM gets an all-zero last row, `n_alpha`. N (-1), and any
	# other index without a PWM row, is sent to it: as uint8, -1 is 255. `X`
	# is always a fresh array here, so the caller's input is not modified.
	n_alpha = motif_pwms.shape[0]
	numpy.minimum(X.view(numpy.uint8), n_alpha, out=X.view(numpy.uint8))

	pwms_n = numpy.zeros((n_alpha + 1, motif_pwms.shape[1]),
		dtype=motif_pwms.dtype, order='F')
	pwms_n[:n_alpha] = motif_pwms

	# Use a fast numba function to run the core algorithm. The hits come back
	# as flat columns ordered by output DataFrame: each motif's forward hits,
	# then, if scanned, the hits of its reverse complement.
	n_ = n_motifs // 2 if reverse_complement else n_motifs
	step = 2 if reverse_complement else 1

	order = numpy.arange(n_motifs, dtype=numpy.int64)
	if reverse_complement:
		order = order.reshape(2, n_).T.flatten()

	# The code of the q letters at every position, shared by every motif.
	q = _qmer_width(n_alpha + 1)
	codes = _qmer_codes(X, n_alpha + 1, q)

	_col_max, _gap, _cuts, _tops = _score_bounds(motif_pwms, motif_lengths,
		_score_thresholds)
	_blk_off, _blk_w, _blk_rest, _blk_start = _block_layout(_col_max, _gap,
		motif_lengths, q)
	# The motif-strands differ in cost, and some are skipped outright, so the
	# scan's prange hands them out one at a time instead of in one contiguous
	# block per thread. The chunk size is set from Python
	# because setting it inside a cached function disables its cache, and it
	# is restored on every exit so that the caller's setting is unchanged.
	previous = numba.set_parallel_chunksize(1)
	try:
		offsets, seqs, starts, ends, scores, pvals = _fast_hits(X, codes, q,
			X_lengths, pwms_n, motif_lengths, _score_thresholds, bin_size,
			_smallest, _score_to_pvals, _score_to_pvals_lengths, _blk_off,
			_blk_w, _blk_rest, _blk_start, _cuts, _tops, order)
	finally:
		numba.set_parallel_chunksize(previous)

	if return_counts == True:
		return numpy.diff(offsets[::step]).astype('int32')

	# Convert the results to pandas DataFrames
	names = ['motif_name', 'motif_idx', 'sequence_name', 'start', 'end', 
		'strand', 'score', 'p-value']

	string_names = motif_names.dtype.kind == 'U'
	if sequence_names is not None:
		sequence_names_ = sequence_names.astype(object)

	# When the names are all strings, every hit goes into one DataFrame, in
	# output order, and each motif's DataFrame is a copy of its rows: one
	# constructor call per fimo call instead of one per motif. Every column
	# gets the values, and the dtype, it gets when each DataFrame is built
	# alone. Any other names are passed per DataFrame as a list, so that pandas
	# infers each motif_name column's dtype from that motif's name alone.
	n_total = offsets[-1]
	if string_names and n_total > 0:
		counts = numpy.diff(offsets[::step])
		labels = numpy.array(list(motif_names[:n_]), dtype=object)
		strands = numpy.array(['+', '-'][:step] * n_, dtype=object)

		if sequence_names is not None:
			sequence_idxs = sequence_names_[seqs]
		else:
			sequence_idxs = seqs

		table = pandas.DataFrame({
			'motif_name': numpy.repeat(labels, counts),
			'motif_idx': numpy.repeat(numpy.arange(n_, dtype=numpy.int64),
				counts),
			'sequence_name': sequence_idxs,
			'start': starts,
			'end': ends,
			'strand': numpy.repeat(strands, numpy.diff(offsets)),
			'score': scores,
			'p-value': pvals
		})

	hits, empty = [], None
	for i in range(n_):
		a, b = offsets[step*i], offsets[step*i + step]
		n_fwd = offsets[step*i + 1] - a
		n = b - a

		# Every DataFrame without hits is the same, so build one, in the way
		# that sets its dtypes, and copy it for the rest. A copy is
		# consolidated, which makes copying it cheaper, so each further one is
		# copied from the last.
		if n == 0:
			if empty is None:
				empty = pandas.DataFrame([], columns=['sequence_name', 
					'start', 'end', 'score', 'p-value'])
				empty['strand'] = []
				empty['motif_name'] = []
				empty['motif_idx'] = numpy.ones(0, dtype='int64')

				if sequence_names is not None:
					empty['sequence_name'] = sequence_names[
						empty['sequence_name']]

				empty = empty[names]
			else:
				empty = empty.copy()

			hits.append(empty)
			continue

		if string_names:
			if n == n_total:
				hits.append(table)
			else:
				hits.append(table.iloc[a:b].reset_index(drop=True))

			continue

		names_ = [motif_names[i]] * n
		strands = numpy.empty(n, dtype=object)
		strands[:n_fwd] = '+'
		strands[n_fwd:] = '-'

		if sequence_names is not None:
			sequence_idxs = sequence_names_[seqs[a:b]]
		else:
			sequence_idxs = seqs[a:b]

		hits.append(pandas.DataFrame({
			'motif_name': names_,
			'motif_idx': numpy.full(n, i, dtype=numpy.int64),
			'sequence_name': sequence_idxs,
			'start': starts[a:b],
			'end': ends[a:b],
			'strand': strands,
			'score': scores[a:b],
			'p-value': pvals[a:b]
		}))

	if dim == 1:
		hits = [df for df in hits if len(df) > 0]
		if len(hits) == 0:
			return []

		hits = pandas.concat(hits)
		_names = numpy.unique(hits['sequence_name'])
		hits = [hits[hits['sequence_name'] == name].reset_index(drop=True)
			for name in _names]

	return hits

