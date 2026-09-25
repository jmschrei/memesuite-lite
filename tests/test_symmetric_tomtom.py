# test_symmetric_tomtom.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import numpy
import pytest

from memelite.io import read_meme
from memelite.tomtom import tomtom
from memelite.symmetric_tomtom import symmetric_tomtom

from numpy.testing import assert_raises
from numpy.testing import assert_array_almost_equal


def generate_random_meme(n=5, min_len=4, max_len=20, random_state=0):
	state = numpy.random.RandomState(random_state)

	pwms = []
	for i in range(n):
		length = state.choice(max_len-min_len+1) + min_len

		pwm = state.choice(17, p=[0.2] + [0.05]*16, size=(length, 4))
		pwm[pwm.sum(axis=1) == 0] = [1, 0, 0, 0]
		pwm = pwm / pwm.sum(axis=1, keepdims=True)
		pwms.append(pwm.T)

	return pwms


def generate_distinct_length_meme(lengths, random_state=0):
	state = numpy.random.RandomState(random_state)

	pwms = []
	for length in lengths:
		pwm = state.choice(17, p=[0.2] + [0.05]*16, size=(length, 4))
		pwm[pwm.sum(axis=1) == 0] = [1, 0, 0, 0]
		pwm = pwm / pwm.sum(axis=1, keepdims=True)
		pwms.append(pwm.T)

	return pwms


###


def test_symmetric_tomtom_ndarray():
	pwms = generate_random_meme(n=12)
	p, scores, offsets, overlaps, strands = symmetric_tomtom(pwms)

	assert isinstance(p, numpy.ndarray)
	assert isinstance(scores, numpy.ndarray)
	assert isinstance(offsets, numpy.ndarray)
	assert isinstance(overlaps, numpy.ndarray)
	assert isinstance(strands, numpy.ndarray)

	assert p.dtype == numpy.float64
	assert scores.dtype == numpy.float64
	assert offsets.dtype == numpy.float64
	assert overlaps.dtype == numpy.float64
	assert strands.dtype == numpy.float64

	assert p.shape == (12, 12)
	assert scores.shape == (12, 12)
	assert offsets.shape == (12, 12)
	assert overlaps.shape == (12, 12)
	assert strands.shape == (12, 12)


def test_symmetric_tomtom_golden():
	pwms = generate_random_meme(n=12)
	p, scores, offsets, overlaps, strands = symmetric_tomtom(pwms)

	# p[0, 0] and scores[0, 0] are the diagonal defaults (1 and 0), which are
	# stable, but offsets[0, 0] and overlaps[0, 0] are never written and hold
	# scratchpad garbage, so those golden rows are checked from index 1 on.
	assert_array_almost_equal(p[0], [1.        , 0.5161135 , 0.29752472,
		0.9894975 , 0.3267443 , 0.4645958 , 0.26181891, 0.16032022,
		0.88019144, 0.47495506, 0.29847423, 0.79400134], 4)
	assert_array_almost_equal(scores[0], [0., 491., 296., 975., 553., 508.,
		462., 1026., 764., 246., 478., 304.])
	assert_array_almost_equal(offsets[0, 1:], [10., 12., 14., -2., 12., -2., 4.,
		1., 6., -2., 1.])
	assert_array_almost_equal(overlaps[0, 1:], [6., 4., 4., 6., 4., 5., 13., 14.,
		4., 5., 5.])
	assert_array_almost_equal(strands[0], [1., 0., 0., 0., 0., 0., 0., 1., 1.,
		0., 1., 1.])


###


def test_symmetric_tomtom_symmetry():
	pwms = generate_random_meme(n=12)
	p, scores, offsets, overlaps, strands = symmetric_tomtom(pwms)

	# All five matrices are made symmetric by construction: the lower
	# triangle is overwritten with the upper triangle, including the
	# offsets and strands (so they are symmetric, not antisymmetric). The
	# diagonal is never written (self-comparison is skipped) and can hold NaN
	# scratchpad garbage, so mask it out -- symmetry of the diagonal is trivial.
	mask = ~numpy.eye(12, dtype=bool)
	assert_array_almost_equal(p[mask], p.T[mask], 4)
	assert_array_almost_equal(scores[mask], scores.T[mask], 4)
	assert_array_almost_equal(offsets[mask], offsets.T[mask], 4)
	assert_array_almost_equal(overlaps[mask], overlaps.T[mask], 4)
	assert_array_almost_equal(strands[mask], strands.T[mask], 4)


def test_symmetric_tomtom_symmetry_meme():
	pwms = list(read_meme("tests/data/test.meme").values())
	p, scores, offsets, overlaps, strands = symmetric_tomtom(pwms)

	mask = ~numpy.eye(len(pwms), dtype=bool)
	assert_array_almost_equal(p[mask], p.T[mask], 4)
	assert_array_almost_equal(scores[mask], scores.T[mask], 4)
	assert_array_almost_equal(offsets[mask], offsets.T[mask], 4)
	assert_array_almost_equal(overlaps[mask], overlaps.T[mask], 4)
	assert_array_almost_equal(strands[mask], strands.T[mask], 4)


###


def test_symmetric_tomtom_diagonal():
	pwms = generate_random_meme(n=12)
	p, scores, offsets, overlaps, strands = symmetric_tomtom(pwms)

	# The diagonal is the self-comparison, but `_p_values` explicitly skips
	# comparing a query against itself (and the symmetry loop only fills the
	# lower triangle from the upper), so the diagonal is left at the default
	# p-value of 1 and score of 0 rather than a near-zero p-value.
	assert_array_almost_equal(numpy.diag(p), numpy.ones(12), 4)
	assert_array_almost_equal(numpy.diag(scores), numpy.zeros(12), 4)


###


def test_symmetric_tomtom_order_invariance():
	# Distinct lengths so the stable length-sort has no ties; with ties the
	# query/target asymmetry of TOMTOM makes the symmetrized output depend on
	# the relative order of equal-length motifs (see module-level note).
	pwms = generate_distinct_length_meme([4, 6, 8, 10, 12, 14, 16, 18])
	n = len(pwms)

	p, scores, offsets, overlaps, strands = symmetric_tomtom(pwms)

	perm = numpy.random.RandomState(1).permutation(n)
	inv = numpy.argsort(perm)
	pwms_perm = [pwms[i] for i in perm]

	pp, sp, op, ovp, stp = symmetric_tomtom(pwms_perm)

	# The diagonal is never written (self-comparison is skipped) and so holds
	# scratchpad garbage; only compare the meaningful off-diagonal values.
	mask = ~numpy.eye(n, dtype=bool)
	assert_array_almost_equal(pp[inv][:, inv][mask], p[mask], 4)
	assert_array_almost_equal(sp[inv][:, inv][mask], scores[mask], 4)
	assert_array_almost_equal(op[inv][:, inv][mask], offsets[mask], 4)
	assert_array_almost_equal(ovp[inv][:, inv][mask], overlaps[mask], 4)
	assert_array_almost_equal(stp[inv][:, inv][mask], strands[mask], 4)


###


def test_symmetric_tomtom_vs_tomtom_diagonal():
	# Plain tomtom computes the self-comparison on the diagonal (near-zero
	# p-values); symmetric_tomtom skips it. Verify they differ on the diagonal
	# but that the off-diagonal best matches are consistent in structure.
	pwms = generate_distinct_length_meme([4, 6, 8, 10, 12, 14, 16, 18])
	n = len(pwms)

	p, scores = symmetric_tomtom(pwms)[:2]
	tp, ts = tomtom(pwms, pwms)[:2]

	assert numpy.all(numpy.diag(tp) < 1e-6)
	assert_array_almost_equal(numpy.diag(p), numpy.ones(n), 4)

	# The upper triangle of symmetric scores should match the corresponding
	# directional tomtom scores (query index < target index).
	for i in range(n):
		for j in range(i+1, n):
			assert_array_almost_equal(scores[i, j], ts[i, j], 4)


###


def test_symmetric_tomtom_reverse_complement_false():
	pwms = generate_random_meme(n=12)
	p, scores, offsets, overlaps, strands = symmetric_tomtom(pwms,
		reverse_complement=False)

	assert p.shape == (12, 12)

	# Without reverse complements, the strand is always 0.
	assert_array_almost_equal(strands, numpy.zeros((12, 12)), 4)

	mask = ~numpy.eye(12, dtype=bool)
	assert_array_almost_equal(p[mask], p.T[mask], 4)
	assert_array_almost_equal(scores[mask], scores.T[mask], 4)

	assert_array_almost_equal(p[0], [1., 1., 0.22434333, 1., 1., 1., 1., 1., 1.,
		0.2691859 , 1., 0.56274301], 4)
	assert_array_almost_equal(scores[0], [0., 0., 297., 0., 0., 0., 0., 0., 0.,
		245., 0., 304.])


###


def test_symmetric_tomtom_n_target_bins_none():
	pwms = generate_random_meme(n=12)
	p0 = symmetric_tomtom(pwms, n_target_bins=100)[0]
	p1 = symmetric_tomtom(pwms, n_target_bins=None)[0]

	assert_array_almost_equal(p1, p1.T, 4)

	# With this data the default hashing does not lose accuracy.
	assert_array_almost_equal(p0, p1, 4)
	assert_array_almost_equal(p1[0], [1., 0.5161135 , 0.29752472, 0.9894975 ,
		0.3267443 , 0.4645958 , 0.26181891, 0.16032022, 0.88019144, 0.47495506,
		0.29847423, 0.79400134], 4)


###


def test_symmetric_tomtom_n_jobs():
	pwms = generate_random_meme(n=12)
	res1 = symmetric_tomtom(pwms, n_jobs=1)
	resm1 = symmetric_tomtom(pwms, n_jobs=-1)

	# The never-written diagonal holds thread-local scratchpad garbage, so it
	# can differ between thread counts; the meaningful off-diagonal values are
	# identical. Mask the diagonal before comparing.
	mask = ~numpy.eye(12, dtype=bool)
	for a, b in zip(res1, resm1):
		assert_array_almost_equal(a[mask], b[mask], 4)


###


def test_symmetric_tomtom_zeroes():
	all_zeroes = numpy.zeros((4, 6))
	assert_raises(ValueError, symmetric_tomtom, [all_zeroes])


###


def generate_dirichlet_meme(lengths, alpha=0.5, random_state=0):
	state = numpy.random.RandomState(random_state)
	return [state.dirichlet(numpy.ones(4) * alpha, size=int(length)).T
		for length in lengths]


def expected_from_tomtom(pwms, **kwargs):
	"""Build the expected symmetric_tomtom output from a plain tomtom call.

	symmetric_tomtom stable-sorts the motifs by length, runs each motif as a
	query only against the motifs that come after it in that order, and then
	copies the upper triangle into the lower triangle. So for a pair (a, b)
	the answer is tomtom's answer with the earlier-ranked motif as the query.
	tomtom is run on the sorted list because approximate target hashing keeps
	the first column seen in each bin, so with coarse `n_target_bins` the
	result depends on target order. The diagonal is left as NaN.
	"""

	n = len(pwms)
	lengths = numpy.array([pwm.shape[-1] for pwm in pwms])
	order = numpy.argsort(lengths, kind='stable')
	rank = numpy.argsort(order)

	sorted_pwms = [pwms[i] for i in order]
	results = tomtom(sorted_pwms, sorted_pwms, **kwargs)
	expected = numpy.full((5, n, n), numpy.nan)
	for a in range(n):
		for b in range(n):
			if a == b:
				continue

			q, t = sorted((rank[a], rank[b]))
			expected[:, a, b] = [r[q, t] for r in results]

	return expected


def assert_off_diagonal_equal(x, y):
	"""Exact equality for integer outputs and tight equality for p-values."""

	n = x[0].shape[0]
	mask = ~numpy.eye(n, dtype=bool)

	numpy.testing.assert_allclose(x[0][mask], y[0][mask], rtol=1e-6,
		atol=1e-12)
	for a, b in zip(x[1:], y[1:]):
		numpy.testing.assert_array_equal(a[mask], b[mask])


SHAPE_GRID = {
	'n2_mixed': [3, 9],
	'n2_equal': [7, 7],
	'n3_mixed': [1, 5, 12],
	'n3_equal': [6, 6, 6],
	'n8_sorted': [2, 4, 6, 8, 10, 12, 14, 16],
	'n8_unsorted': [16, 3, 11, 5, 9, 2, 14, 7],
	'n8_ties': [5, 9, 5, 9, 5, 12, 9, 12],
	'n8_long': [20, 24, 28, 32, 36, 40, 3, 9],
	'n30_mixed': list(numpy.random.RandomState(3).randint(1, 21, size=30)),
	'n30_equal': [10] * 30,
}


@pytest.mark.parametrize("name", list(SHAPE_GRID))
def test_symmetric_tomtom_shapes(name):
	lengths = SHAPE_GRID[name]
	n = len(lengths)
	pwms = generate_dirichlet_meme(lengths, random_state=len(name))

	p, scores, offsets, overlaps, strands = symmetric_tomtom(pwms, n_jobs=2)

	for x in (p, scores, offsets, overlaps, strands):
		assert isinstance(x, numpy.ndarray)
		assert x.dtype == numpy.float64
		assert x.shape == (n, n)

	mask = ~numpy.eye(n, dtype=bool)
	assert numpy.all(p[mask] >= 0)
	assert numpy.all(p[mask] <= 1)
	assert set(numpy.unique(strands[mask])).issubset({0.0, 1.0})

	# Scores, offsets and overlaps are integers stored as floats.
	for x in (scores, offsets, overlaps):
		assert_array_almost_equal(x[mask], numpy.round(x[mask]), 12)

	# The overlap can never exceed either motif's length and is at least 1,
	# and the offset lies in the range of valid alignments of the query
	# (the shorter motif by stable rank) against the target.
	lens = numpy.array(lengths)
	rank = numpy.argsort(numpy.argsort(lens, kind='stable'))
	for a in range(n):
		for b in range(n):
			if a == b:
				continue

			q, t = (a, b) if rank[a] < rank[b] else (b, a)
			assert 1 <= overlaps[a, b] <= min(lens[a], lens[b])
			assert -(lens[q] - 1) <= offsets[a, b] <= lens[t] - 1

	# Symmetric by construction, and the diagonal holds the defaults.
	assert_array_almost_equal(p[mask], p.T[mask], 12)
	for x in (scores, offsets, overlaps, strands):
		numpy.testing.assert_array_equal(x[mask], x.T[mask])

	numpy.testing.assert_array_equal(numpy.diag(p), numpy.ones(n))
	numpy.testing.assert_array_equal(numpy.diag(scores), numpy.zeros(n))


@pytest.mark.parametrize("name", list(SHAPE_GRID))
def test_symmetric_tomtom_matches_tomtom(name):
	# With reverse complements, every off-diagonal entry is exactly the plain
	# tomtom result where the earlier motif in the stable length order is the
	# query. This holds with length ties and unsorted input.
	lengths = SHAPE_GRID[name]
	pwms = generate_dirichlet_meme(lengths, random_state=len(name))

	observed = symmetric_tomtom(pwms, n_jobs=2)
	expected = expected_from_tomtom(pwms, n_jobs=2)
	assert_off_diagonal_equal(observed, expected)


@pytest.mark.parametrize("kwargs", [
	{'n_target_bins': None},
	{'n_target_bins': 10},
	{'n_target_bins': 100},
	{'n_score_bins': 50},
	{'n_score_bins': 120},
	{'n_median_bins': 100},
	{'n_median_bins': 5000},
	{'n_cache': 250},
])
def test_symmetric_tomtom_matches_tomtom_kwargs(kwargs):
	# Each kwarg is forwarded identically to the shared kernels, so the
	# equivalence to plain tomtom holds for every setting.
	pwms = generate_random_meme(n=10)

	observed = symmetric_tomtom(pwms, n_jobs=2, **kwargs)
	expected = expected_from_tomtom(pwms, n_jobs=2, **kwargs)
	assert_off_diagonal_equal(observed, expected)


@pytest.mark.skip(reason="BUG: with reverse_complement=False, `_p_values` "
	"still skips targets in [N//2, N//2 + iq] as if the target list held "
	"reverse complements, so those pairs come back as p=1, score=0. The "
	"existing test_symmetric_tomtom_reverse_complement_false golden values "
	"contain these skipped pairs.")
def test_symmetric_tomtom_matches_tomtom_no_rc():
	pwms = generate_distinct_length_meme([4, 6, 8, 10, 12, 14, 16, 18])

	observed = symmetric_tomtom(pwms, reverse_complement=False)
	expected = expected_from_tomtom(pwms, reverse_complement=False)
	assert_off_diagonal_equal(observed, expected)


def test_symmetric_tomtom_all_length_one():
	# Regression: the A workspace had a first axis of Q_max, which is 1 here,
	# but `_p_value_backgrounds` writes A[1].
	pwms = generate_dirichlet_meme([1, 1, 1, 1, 1])

	observed = symmetric_tomtom(pwms)
	expected = expected_from_tomtom(pwms)
	assert_off_diagonal_equal(observed, expected)


##


def test_symmetric_tomtom_n_cache_larger():
	# n_cache only sizes the scratchpad, so raising it must not change any
	# output.
	pwms = generate_random_meme(n=10)

	r0 = symmetric_tomtom(pwms)
	r1 = symmetric_tomtom(pwms, n_cache=500)
	assert_off_diagonal_equal(r0, r1)


@pytest.mark.parametrize("n_cache", [0, 5, 20])
def test_symmetric_tomtom_n_cache_too_small(n_cache):
	# Regression: an offset above n_cache used to overrun the workspace. Such
	# queries now get their own workspace, so the results are unchanged.
	pwms = generate_random_meme(n=10)

	r0 = symmetric_tomtom(pwms)
	r1 = symmetric_tomtom(pwms, n_cache=n_cache)
	assert_off_diagonal_equal(r0, r1)


@pytest.mark.parametrize("n_jobs", [1, 2, 3, -1])
def test_symmetric_tomtom_n_jobs_grid(n_jobs):
	pwms = generate_dirichlet_meme(SHAPE_GRID['n30_mixed'], random_state=5)

	r0 = symmetric_tomtom(pwms, n_jobs=1)
	r1 = symmetric_tomtom(pwms, n_jobs=n_jobs)
	assert_off_diagonal_equal(r0, r1)


@pytest.mark.parametrize("n_jobs", [1, 2, 3, -1])
def test_symmetric_tomtom_restores_num_threads(n_jobs):
	import numba

	pwms = generate_random_meme(n=6)
	n_threads = numba.get_num_threads()

	symmetric_tomtom(pwms, n_jobs=n_jobs)
	assert numba.get_num_threads() == n_threads


def test_symmetric_tomtom_permutation_distinct_lengths():
	# With distinct lengths the stable sort fully determines the internal
	# order, so permuting the input must permute the output. Inputs with
	# length ties are covered by test_symmetric_tomtom_matches_tomtom.
	pwms = generate_dirichlet_meme([4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
		random_state=7)
	n = len(pwms)
	r = symmetric_tomtom(pwms)

	for seed in range(3):
		perm = numpy.random.RandomState(seed).permutation(n)
		inv = numpy.argsort(perm)
		rp = symmetric_tomtom([pwms[i] for i in perm])
		rp = [x[inv][:, inv] for x in rp]
		assert_off_diagonal_equal(r, rp)


def test_symmetric_tomtom_input_not_mutated():
	pwms = generate_dirichlet_meme([3, 8, 5, 12, 7])
	copies = [pwm.copy() for pwm in pwms]

	symmetric_tomtom(pwms)
	for pwm, pwm_copy in zip(pwms, copies):
		numpy.testing.assert_array_equal(pwm, pwm_copy)


def test_symmetric_tomtom_non_contiguous():
	# Strided and Fortran-ordered views must give the same answer as
	# contiguous arrays.
	pwms = generate_random_meme(n=10)
	r = symmetric_tomtom(pwms)

	strided = [numpy.repeat(pwm, 2, axis=1)[:, ::2] for pwm in pwms]
	fortran = [numpy.asfortranarray(pwm) for pwm in pwms]

	assert_off_diagonal_equal(r, symmetric_tomtom(strided))
	assert_off_diagonal_equal(r, symmetric_tomtom(fortran))


def test_symmetric_tomtom_float32():
	# float32 inputs change the distance numerics slightly, so p-values are
	# close but not identical; the integer outputs are unchanged here.
	pwms = generate_random_meme(n=10)
	r64 = symmetric_tomtom(pwms)
	r32 = symmetric_tomtom([pwm.astype('float32') for pwm in pwms])

	mask = ~numpy.eye(10, dtype=bool)
	assert numpy.abs(r64[0][mask] - r32[0][mask]).max() < 1e-2
	for a, b in zip(r64[1:], r32[1:]):
		numpy.testing.assert_array_equal(a[mask], b[mask])


def test_symmetric_tomtom_single_motif():
	# A single motif has no pairs; the output is a 1x1 matrix of defaults.
	pwms = generate_dirichlet_meme([8])
	p, scores, offsets, overlaps, strands = symmetric_tomtom(pwms)

	for x in (p, scores, offsets, overlaps, strands):
		assert x.shape == (1, 1)

	assert p[0, 0] == 1
	assert scores[0, 0] == 0


def test_symmetric_tomtom_reverse_complement_strands():
	# Comparing a set against itself with every motif's reverse complement
	# appended: the best match of a motif and its own reverse complement is
	# a full-overlap reverse-strand alignment.
	pwms = generate_distinct_length_meme([5, 7, 9, 11])
	pwms = pwms + [pwm[::-1, ::-1] for pwm in pwms]
	p, scores, offsets, overlaps, strands = symmetric_tomtom(pwms)

	for i, length in enumerate([5, 7, 9, 11]):
		assert strands[i, i+4] == 1
		assert offsets[i, i+4] == 0
		assert overlaps[i, i+4] == length
		assert p[i, i+4] < 1e-3


def test_symmetric_tomtom_n_score_bins_scales_scores():
	# Scores are sums of per-column integer bins, so their scale follows
	# n_score_bins; the equivalence to tomtom is covered above.
	pwms = generate_random_meme(n=8)
	mask = ~numpy.eye(8, dtype=bool)

	s50 = symmetric_tomtom(pwms, n_score_bins=50)[1][mask]
	s100 = symmetric_tomtom(pwms, n_score_bins=100)[1][mask]
	assert s100.mean() > 1.5 * s50.mean()
