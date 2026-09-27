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

	n, l = log_pwm.shape

	log_bg = math.log2(0.25)
	int_log_pwm = numpy.round(log_pwm / bin_size).astype(numpy.int32)

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

	logpdf = numpy.empty(largest - smallest + 1)
	old_logpdf = -numpy.inf * numpy.ones(largest - smallest + 1)
	for i in range(n):
		idx = int_log_pwm[i, 0] - smallest
		old_logpdf[idx] = logaddexp2(old_logpdf[idx], log_bg)

	for i in range(1, l):
		for j in range(largest - smallest + 1):
			logpdf[j] = -numpy.inf

		for j, x in enumerate(old_logpdf):
			if x != -numpy.inf:
				for k in range(n):
					idx = j + int_log_pwm[k, i]
					logpdf[idx] = logaddexp2(logpdf[idx], log_bg + x)

		for j in range(largest - smallest + 1):
			old_logpdf[j] = logpdf[j]

	# `old_logpdf` holds the distribution over every column. `logpdf` is only
	# written by the loop above, which does not run for a single-column PWM.
	logpdf = old_logpdf

	for i in range(len(logpdf) - 2, -1, -1):
		logpdf[i] = logaddexp2(logpdf[i], logpdf[i + 1])

	return smallest, logpdf


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
# weights of q consecutive columns for every code, so q columns cost one
# lookup. q is the largest value up to `_QMAX` whose tables have at most
# `_TABLE_MAX` entries: 5 for DNA, whose 4 letters plus N give 5**5 = 3125.
# A window is tested once after the first two groups of q columns, then after
# every group. Chosen by measurement on the benchmark: two groups of 5 beat
# two or three groups of 4 or 3, and a test after every column was slower
# than no test at all, because the branch that leaves the window mispredicts.
# Codes are uint16, so `_TABLE_MAX` must stay at most 65536.
_QMAX = 5
_TABLE_MAX = 3125


def _score_bounds(pwm, pwm_lengths, thresholds):
	"""Upper bounds that let `_fast_hits` abandon a window early.

	A column can add at most its largest entry to a window's score, or 0.0 when
	the position is an N, which reads the kernel's all-zero row. For global
	column c = s + j of the motif whose columns start at s, `rest[c + 1]` is
	the most that columns j+1 onwards can still add, so `rest[s + p]` bounds
	what remains after the
	first p columns (`rest[s + n]` and `rest[0]` are 0.0). `tops[k]` is the
	largest score motif k can reach at all.

	A window is abandoned only when its partial score plus the remaining bound
	is at most `cuts[k]`, the threshold minus a margin. The partial score is a
	sum of table entries that each sum q columns (with +0.0 for columns past
	the motif's end, which is exact), so it is one more order of addition of
	the same terms. Every float sum
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
		col_max = pwm.max(axis=0)
		col_max = numpy.where(numpy.isnan(col_max), numpy.inf,
			numpy.maximum(col_max, 0.0))
		col_abs = numpy.where(numpy.isfinite(pwm), numpy.abs(pwm), 0.0).max(
			axis=0)

	# One row per motif, zero-padded, so each suffix sums only its own motif.
	padded = numpy.zeros((n_motifs, max_width + 1))
	padded[:, :max_width][mask] = col_max
	suffix = numpy.cumsum(padded[:, ::-1], axis=1)[:, ::-1]

	rest = numpy.zeros(int(lengths[-1]) + 1)
	rest[1:] = suffix[:, 1:][mask]
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
	return rest, cuts, tops


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


@numba.njit(inline='always')
def _qmer_tables(tab, pwm, off, n, q, stride, n_groups):
	"""tab[g * stride + c] = the summed weights of columns g*q .. g*q + q-1
	of the motif whose columns start at `off`, for the q letters with code c.
	`stride` is at least the number of codes, n_rows**q. Columns at or past
	`n` weigh 0.0. Each group's table is built in place, from its last column
	to its first, each column multiplying its size by the number of rows."""

	n_rows = numpy.int64(pwm.shape[0])
	q = numpy.int64(q)
	n = numpy.int64(n)
	off = numpy.int64(off)
	stride = numpy.int64(stride)
	for g in range(numpy.int64(n_groups)):
		t0 = g * stride
		tab[t0] = 0.0
		size = numpy.int64(1)
		for d in range(q - 1, -1, -1):
			col = g * q + d
			for t in range(size - 1, -1, -1):
				v = tab[t0 + t]
				for s in range(n_rows - 1, -1, -1):
					w = 0.0
					if col < n:
						w = numpy.float64(pwm[s, off + col])
					tab[t0 + t * n_rows + s] = v + w
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
	bin_size, smallest, score_to_pvals, score_to_pval_lengths, rest, cuts, tops,
	order):
	"""Scan every motif over every sequence and return the hits as columns.

	`codes` and `q` come from `_qmer_codes`: the bound that abandons a window
	reads q columns per table lookup. The hits of each motif are collected in
	a list while it is scanned, then copied into flat arrays in the order given
	by `order`, so that the hits of motif `order[t]` are rows `offsets[t]` to
	`offsets[t+1]` of the returned sequence index, start, end, score and
	p-value arrays.
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
		k = numpy.uint64(k)
		thresh = score_threshold[k]
		cut = cuts[k]
		off = numpy.uint64(pwm_lengths[k])
		p = numpy.uint64(min(numpy.uint64(2) * q, n))

		# Every group's table starts `_TABLE_MAX` entries after the last,
		# whatever the number of codes, so the second lookup of the prefix is a
		# constant offset from the table's start and needs no register of its
		# own. It is set inside the prange body because a value set outside is
		# passed into the body as an argument, which LLVM cannot fold.
		stride = numpy.uint64(_TABLE_MAX)

		# Skip a motif whose best possible score cannot pass its threshold.
		if tops[k] > cut:
			# One table per group of q columns, and at least two, so that the
			# first test always reads two. Groups past the motif's end are
			# all 0.0. An N reads the all-zero last row of `pwm` and adds
			# +0.0, in the tables and in the rescoring loop below.
			n_groups = max((n + q - numpy.uint64(1)) // q, numpy.uint64(2))
			tab = numpy.empty(n_groups * stride, dtype=numpy.float64)
			_qmer_tables(tab, pwm, off, n, q, stride, n_groups)
			rest_p = rest[off+p]

			for l in range(n_chroms):        
				start = numpy.uint64(chrom_lengths[l])
				end = numpy.uint64(chrom_lengths[l+1])
				
				# Windows start at positions x from `start` to `xe` - 1.
				m = numpy.uint64(max(numpy.int64(end - start) -
					numpy.int64(n) + 1, 0))
				xe = start + m
				x = start
				while x < xe:
					# The first 2q columns in two lookups, in a loop of its own
					# that moves to the next window until one survives this
					# test, which about 1% of windows do on real motifs. With
					# nothing else in the loop, its pointers, position and
					# bounds stay in registers; in one loop with the code
					# below, LLVM reloaded and spilled them for every window.
					# Only the bound test uses these columns.
					bound = 0.0
					while x < xe:
						bound = tab[numpy.uint64(codes[x])] + tab[stride +
							numpy.uint64(codes[x+q])]
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
					while g < n_groups:
						bound += tab[g * stride + numpy.uint64(codes[base+g*q])]
						j1 = min((g + numpy.uint64(1)) * q, n)
						if bound + rest[off+j1] <= cut:
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

	_smallest, _score_to_pvals = _all_pwm_to_mapping(motif_pwms, motif_lengths, 
		bin_size)
	_score_to_pvals_lengths = [0]
	_score_thresholds = numpy.empty(n_motifs, dtype=numpy.float32)

	for i in range(n_motifs):	
		_score_to_pvals_lengths.append(len(_score_to_pvals[i]))

		idx = numpy.where(_score_to_pvals[i] < log_threshold)[0]
		if len(idx) > 0:
			_score_thresholds[i] = (idx[0] + _smallest[i]) * bin_size                              
		else:
			_score_thresholds[i] = float("inf")

	_score_to_pvals = numpy.concatenate(_score_to_pvals)
	_score_to_pvals_lengths = numpy.cumsum(_score_to_pvals_lengths)

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
		X = ((sequences.argmax(axis=1) + 1) * sequences.sum(axis=1)) - 1
		X_lengths = numpy.arange(X.shape[0]+1) * X.shape[-1]

		X = X.astype(numpy.int8).flatten()
		X_lengths = X_lengths.astype(numpy.int64)

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

	_rest, _cuts, _tops = _score_bounds(motif_pwms, motif_lengths,
		_score_thresholds)
	offsets, seqs, starts, ends, scores, pvals = _fast_hits(X, codes, q,
		X_lengths, pwms_n, motif_lengths, _score_thresholds, bin_size,
		_smallest, _score_to_pvals, _score_to_pvals_lengths, _rest, _cuts,
		_tops, order)

	if return_counts == True:
		return numpy.diff(offsets[::step]).astype('int32')

	# Convert the results to pandas DataFrames
	names = ['motif_name', 'motif_idx', 'sequence_name', 'start', 'end', 
		'strand', 'score', 'p-value']

	# Names that are all strings go in an object column directly; anything
	# else is passed as a list so pandas infers the column's dtype.
	string_names = motif_names.dtype.kind == 'U'
	if sequence_names is not None:
		sequence_names_ = sequence_names.astype(object)

	hits, empty = [], None
	for i in range(n_):
		a, b = offsets[step*i], offsets[step*i + step]
		n_fwd = offsets[step*i + 1] - a
		n = b - a

		# Every DataFrame without hits is the same, so build one, in the way
		# that sets its dtypes, and copy it for the rest.
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
				hits.append(empty)
			else:
				hits.append(empty.copy())

			continue

		if string_names:
			names_ = numpy.empty(n, dtype=object)
			names_[:] = motif_names[i]
		else:
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

