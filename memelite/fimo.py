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


# How often `_fast_hits` tests whether a window can still reach its threshold:
# once after the first `_PREFIX` columns, then after every `_STEP` columns.
# Chosen by measurement on the benchmark; a test after every column was slower
# than no test at all, because the branch that leaves the window mispredicts.
_PREFIX = 8
_STEP = 4


def _score_bounds(pwm, pwm_lengths, thresholds):
	"""Upper bounds that let `_fast_hits` abandon a window early.

	A column can add at most its largest entry to a window's score, or 0.0 when
	the position is an N, which adds nothing. For global column c = s + j of
	the motif whose columns start at s, `rest[c + 1]` is the most that columns
	j+1 onwards can still add, so `rest[s + p]` bounds what remains after the
	first p columns (`rest[s + n]` and `rest[0]` are 0.0). `tops[k]` is the
	largest score motif k can reach at all.

	A window is abandoned only when its partial score plus the remaining bound
	is at most `cuts[k]`, the threshold minus a margin. Every float sum
	involved (the partial score, the bounds, and the full score in whatever
	order fastmath adds it) has at most n terms of magnitude at most
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
def _fast_hits(X, chrom_lengths, pwm, pwm_lengths, score_threshold, bin_size, 
	smallest, score_to_pvals, score_to_pval_lengths, rest, cuts, tops, order):
	"""Scan every motif over every sequence and return the hits as columns.

	The hits of each motif are collected in a list while it is scanned, then
	copied into flat arrays in the order given by `order`, so that the hits of
	motif `order[t]` are rows `offsets[t]` to `offsets[t+1]` of the returned
	sequence index, start, end, score and p-value arrays.
	"""

	n_motifs = len(pwm_lengths) - 1
	n_chroms = len(chrom_lengths) - 1

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
		p = numpy.uint64(min(numpy.uint64(_PREFIX), n))

		# Skip a motif whose best possible score cannot pass its threshold.
		if tops[k] > cut:
			for l in range(n_chroms):        
				start = numpy.uint64(chrom_lengths[l])
				end = numpy.uint64(chrom_lengths[l+1])
				
				for i in range(end-start-n+1):
					i = numpy.uint64(i)
					base = start + i

					# The first p columns as two interleaved partial sums, which
					# halves the chain of dependent additions. Only the bound
					# test uses them.
					a = 0.0
					b = 0.0
					j = numpy.uint64(0)
					while j + numpy.uint64(1) < p:
						idx = X[base+j]
						if idx != -1:
							a += pwm[numpy.uint64(idx), off+j]
						idx = X[base+j+numpy.uint64(1)]
						if idx != -1:
							b += pwm[numpy.uint64(idx), off+j+numpy.uint64(1)]
						j += numpy.uint64(2)
					if j < p:
						idx = X[base+j]
						if idx != -1:
							a += pwm[numpy.uint64(idx), off+j]

					bound = a + b
					if bound + rest[off+p] <= cut:
						continue

					alive = True
					j0 = p
					while j0 < n:
						j1 = min(j0 + numpy.uint64(_STEP), n)
						for j in range(j0, j1):
							j = numpy.uint64(j)
							idx = X[base+j]
							if idx != -1:
								bound += pwm[numpy.uint64(idx), off+j]

						if bound + rest[off+j1] <= cut:
							alive = False
							break
						j0 = j1

					if not alive:
						continue

					# A window that survives is scored exactly as before, so its
					# score and the decision below are unchanged.
					score = 0.0
					for j in range(n):
						j = numpy.uint64(j)
						
						idx = X[start+i+j]
						if idx == -1:
							continue

						m_idx = numpy.uint64(j + pwm_lengths[k])
						idx = numpy.uint64(idx)
						score += pwm[idx, m_idx]

					if score > thresh:
						score_idx = int(score / bin_size) - smallest[k]                    
						score_idx += score_to_pval_lengths[k]
						hits[k].append((numpy.int64(l), i, i+n, score, 
							2.0 ** score_to_pvals[score_idx]))

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

	# Use a fast numba function to run the core algorithm. The hits come back
	# as flat columns ordered by output DataFrame: each motif's forward hits,
	# then, if scanned, the hits of its reverse complement.
	n_ = n_motifs // 2 if reverse_complement else n_motifs
	step = 2 if reverse_complement else 1

	order = numpy.arange(n_motifs, dtype=numpy.int64)
	if reverse_complement:
		order = order.reshape(2, n_).T.flatten()

	_rest, _cuts, _tops = _score_bounds(motif_pwms, motif_lengths,
		_score_thresholds)
	offsets, seqs, starts, ends, scores, pvals = _fast_hits(X, X_lengths, 
		motif_pwms, motif_lengths, _score_thresholds, bin_size, _smallest, 
		_score_to_pvals, _score_to_pvals_lengths, _rest, _cuts, _tops, order)

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

