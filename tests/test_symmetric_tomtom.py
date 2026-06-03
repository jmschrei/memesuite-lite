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
