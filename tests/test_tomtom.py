# test_tomtom.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import numba
import numpy
import pytest
import pandas

from memelite.io import read_meme
from memelite.tomtom import _binned_median
from memelite.tomtom import _binned_median_z
from memelite.tomtom import _binned_median_block4
from memelite.tomtom import _pairwise_max
from memelite.tomtom import _merge_rc_results
from memelite.tomtom import _p_value_backgrounds
from memelite.tomtom import _p_values
from memelite.tomtom import tomtom

from ._golden_inputs import one_hot_pwms

from numpy.testing import assert_raises
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal

# The grids use two threads to exercise the parallel path while keeping the
# per-thread scratchpad small; n_jobs above numba's thread count is an error.
_N_JOBS = min(2, numba.config.NUMBA_NUM_THREADS)


def _require_threads(n_jobs):
	if n_jobs > numba.config.NUMBA_NUM_THREADS:
		pytest.skip("needs {} numba threads".format(n_jobs))



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


def random_one_hot(shape, probs=None, dtype='int8', random_state=None):
	if not isinstance(shape, tuple) or len(shape) != 3:
		raise ValueError("Shape must be a tuple with 3 dimensions.")

	if not isinstance(random_state, numpy.random.RandomState):
		random_state = numpy.random.RandomState(random_state)

	if isinstance(probs, list):
		probs = numpy.array(probs)
		
	n = shape[1]
	ohe = numpy.zeros(shape, dtype=dtype)

	for i in range(ohe.shape[0]):
		if probs is None:
			probs_ = None
		elif probs.ndim == 1:
			probs_ = probs
		elif probs.shape[0] == 1:
			probs_ = probs[0]
		else:
			probs_ = probs[i] 

		choices = random_state.choice(n, size=shape[2], p=probs_)
		ohe[i, choices, numpy.arange(shape[2])] = 1 

	return one


###


def test_binned_median_odd_short():
	X = numpy.array([0, 4, 2, 1, 3], dtype='float64')
	bins = numpy.zeros((5, 2), dtype='float64')
	counts = numpy.ones(5, dtype='int64')
	assert _binned_median(X, bins, 0, 4, counts) == 2


def test_binned_median_odd_counts():
	X = numpy.array([0, 4, 2, 1, 3], dtype='float64')
	bins = numpy.zeros((5, 2), dtype='float64')
	counts = numpy.ones(5, dtype='int64')
	counts[1] = 6
	assert _binned_median(X, bins, 0, 4, counts) == 4


def test_binned_median_odd_long():
	X = numpy.array([0, 2, 8, 3, 4, 7, 9, 10, 4, 2, 1], dtype='float64')
	bins = numpy.zeros((len(X), 2), dtype='float64')
	counts = numpy.ones(len(X), dtype='int64')
	assert _binned_median(X, bins, 0, 10, counts) == 4


def test_binned_median_even():
	X = numpy.array([0, 1, 2, 3], dtype='float64')
	bins = numpy.zeros((4, 2), dtype='float64')
	counts = numpy.ones(4, dtype='int64')
	assert _binned_median(X, bins, 0, 3, counts) == 1


def test_binned_median_even_counts():
	X = numpy.array([0, 1, 2, 3], dtype='float64')
	bins = numpy.zeros((4, 2), dtype='float64')
	counts = numpy.ones(4, dtype='int64')
	counts[2] = 2
	assert _binned_median(X, bins, 0, 3, counts) == 2


def test_binned_median_long():
	X = numpy.random.RandomState(0).randn(200000)

	bins = numpy.zeros((1000, 2), dtype='float64')
	counts = numpy.ones(len(X), dtype='int64')
	assert_array_almost_equal([_binned_median(X, bins, X.min(), X.max(), 
		counts)], [numpy.median(X)], 2)


###


def test_pairwise_max():
	x = numpy.abs(numpy.random.RandomState(0).randn(100))
	x = x / x.sum()

	y = numpy.abs(numpy.random.RandomState(0).randn(100))
	y = y / y.sum()

	y_csum = numpy.cumsum(y, axis=-1)
	x_csum = numpy.cumsum(x, axis=-1)

	z = numpy.empty(100)
	_pairwise_max(x, y, y_csum, z, 100)

	assert_array_almost_equal(z[:10], [0.000475, 0.00024 , 0.000792, 0.002914, 
		0.003599, 0.002307, 0.002523, 0.000427, 0.000295, 0.001207])
	assert_array_almost_equal(z, x * y_csum + y * x_csum - x * y)


def test_pairwise_max_fallback():
	x = numpy.abs(numpy.random.RandomState(0).randn(100))
	x = x / x.sum()
	x[0] = -1

	y = numpy.abs(numpy.random.RandomState(0).randn(100))
	y = y / y.sum()

	y_csum = numpy.cumsum(y, axis=-1)

	z = numpy.empty(100)
	_pairwise_max(x, y, y_csum, z, 100)

	assert_array_almost_equal(z, y)


###


def test_merge_rc_results():
	results = numpy.random.RandomState(0).randn(1000, 5)

	idxs = results[:500, 1] > results[500:, 1]
	best_p = 1 - (1 - numpy.minimum(results[:500, 0], results[500:, 0])) ** 2
	best_scores = numpy.where(idxs, results[:500, 1], results[500:, 1])
	best_offsets = numpy.where(idxs, results[:500, 2], results[500:, 2])
	best_overlaps = numpy.where(idxs, results[:500, 3], results[500:, 3])

	_merge_rc_results(results)

	assert_array_almost_equal(results[:500, 0], best_p, 4)
	assert_array_almost_equal(results[:500, 1], best_scores, 4)
	assert_array_almost_equal(results[:500, 2], best_offsets, 4)
	assert_array_almost_equal(results[:500, 3], best_overlaps, 4)
	assert_array_almost_equal((1 - results[:500, 4]).astype(bool), idxs)


###


def test_tomtom_ndarray():
	pwms = generate_random_meme(n=20)
	p, scores, offsets, overlaps, strands = tomtom(pwms, pwms)

	assert isinstance(p, numpy.ndarray)
	assert isinstance(scores, numpy.ndarray)
	assert isinstance(offsets, numpy.ndarray)
	assert isinstance(overlaps, numpy.ndarray)
	assert isinstance(strands, numpy.ndarray)

	assert p.dtype == numpy.float64
	assert scores.dtype == numpy.float64
	assert offsets.dtype == numpy.float64
	assert overlaps.dtype == numpy.float64
	assert overlaps.dtype == numpy.float64

	assert p.shape == (20, 20)
	assert scores.shape == (20, 20)
	assert offsets.shape == (20, 20)
	assert overlaps.shape == (20, 20)
	assert strands.shape == (20, 20)

	assert_array_almost_equal(p[0], [1.95399252e-14, 4.04372666e-01, 
		3.33899770e-01, 9.75969973e-01, 7.19287721e-01, 2.02698438e-01, 
		6.09413729e-01, 1.05202910e-01, 9.97854019e-01, 4.39470284e-01, 
		6.39529870e-01, 4.01256302e-01, 9.26477704e-01, 6.25608933e-01, 
		7.54828358e-01, 7.37538119e-01, 9.49853411e-01, 4.77573695e-01, 
		9.66135125e-01, 7.13675581e-01], 6)
	assert_array_almost_equal(scores[0], [1399., 1002.,  991.,  993.,  994., 
		1011.,  995., 1047.,  979.,  988.,  994.,  994., 986., 1015., 1003.,  
		998.,  982., 1023.,  976.,  981.])
	assert_array_almost_equal(offsets[0], [ 0., -4., -6., 14.,  2., -6.,  2.,  
		1., 11., -6., -7., -7.,  5., -4., -6.,  2.,  3.,  1., -7., -8.])
	assert_array_almost_equal(overlaps[0], [16.,  7.,  4.,  4.,  6.,  7.,  5., 
		16.,  3.,  4.,  7.,  5.,  4., 12., 10.,  8.,  5., 16., 6.,  4.])
	assert_array_almost_equal(strands[0], [0., 0., 0., 0., 0., 0., 0., 1., 0., 
		0., 0., 0., 1., 0., 1., 0., 1., 0., 1., 1])


'''
def test_tomtom_pytorch():
	pwms = generate_random_meme(n=20)
	pwms = [torch.from_numpy(pwm) for pwm in pwms]
	p, scores, offsets, overlaps, strands = tomtom(pwms, pwms)

	assert isinstance(p, numpy.ndarray)
	assert isinstance(scores, numpy.ndarray)
	assert isinstance(offsets, numpy.ndarray)
	assert isinstance(overlaps, numpy.ndarray)
	assert isinstance(strands, numpy.ndarray)

	assert p.dtype == numpy.float64
	assert scores.dtype == numpy.float64
	assert offsets.dtype == numpy.float64
	assert overlaps.dtype == numpy.float64
	assert overlaps.dtype == numpy.float64

	assert p.shape == (20, 20)
	assert scores.shape == (20, 20)
	assert offsets.shape == (20, 20)
	assert overlaps.shape == (20, 20)
	assert strands.shape == (20, 20)

	assert_array_almost_equal(p[0], [1.95399252e-14, 4.04372666e-01, 
		3.33899770e-01, 9.75969973e-01, 7.19287721e-01, 2.02698438e-01, 
		6.09413729e-01, 1.05202910e-01, 9.97854019e-01, 4.39470284e-01, 
		6.39529870e-01, 4.01256302e-01, 9.26477704e-01, 6.25608933e-01, 
		7.54828358e-01, 7.37538119e-01, 9.49853411e-01, 4.77573695e-01, 
		9.66135125e-01, 7.13675581e-01], 6)
	assert_array_almost_equal(scores[0], [1399., 1002.,  991.,  993.,  994., 
		1011.,  995., 1047.,  979.,  988.,  994.,  994., 986., 1015., 1003.,  
		998.,  982., 1023.,  976.,  981.])
	assert_array_almost_equal(offsets[0], [ 0., -4., -6., 14.,  2., -6.,  2.,  
		1., 11., -6., -7., -7.,  5., -4., -6.,  2.,  3.,  1., -7., -8.])
	assert_array_almost_equal(overlaps[0], [16.,  7.,  4.,  4.,  6.,  7.,  5., 
		16.,  3.,  4.,  7.,  5.,  4., 12., 10.,  8.,  5., 16., 6.,  4.])
	assert_array_almost_equal(strands[0], [0., 0., 0., 0., 0., 0., 0., 1., 0., 
		0., 0., 0., 1., 0., 1., 0., 1., 0., 1., 1])
'''


@pytest.mark.xfail(strict=True, reason="When the forward and reverse "
	"strand best scores tie, `_merge_rc_results` resolves to strand 1 (`<=`), "
	"so reversing the targets changes the reported strand and offset for tied "
	"pairs. See test_tomtom_reverse_complement_targets_untied.")
def test_tomtom_reverse_complement_targets():
	pwms = generate_random_meme(n=20)
	p0, scores0, offsets0, overlaps0, strands0 = tomtom(pwms, pwms)
	p1, scores1, offsets1, overlaps1, strands1 = tomtom(pwms,
		[p[::-1, ::-1] for p in pwms])

	assert_array_almost_equal(p0, p1, 4)
	assert_array_almost_equal(scores0, scores1, 4)
	assert_array_almost_equal(offsets0, offsets1, 4)
	assert_array_almost_equal(overlaps0, overlaps1, 4)
	assert_array_almost_equal(strands0, 1-strands1)


def test_tomtom_zeroes():
	pwms = generate_random_meme(n=5)
	all_zeroes = numpy.array([
		[0, 0, 0, 0],
		[0, 0, 0, 0],
		[0, 0, 0, 0],
		[0, 0, 0, 0],
		[0, 0, 0, 0]
	])

	assert_raises(ValueError, tomtom, [all_zeroes], pwms)
	assert_raises(ValueError, tomtom, pwms, [all_zeroes])


def test_tomtom_subsets():
	pwms = generate_random_meme(n=20)
	p, scores, offsets, overlaps, strands = tomtom(pwms[:2], pwms)

	assert p.shape == (2, 20)
	assert scores.shape == (2, 20)
	assert offsets.shape == (2, 20)
	assert overlaps.shape == (2, 20)
	assert strands.shape == (2, 20)

	p2, scores2, offsets2, overlaps2, strands2 = tomtom(pwms[:5], pwms)

	assert p2.shape == (5, 20)
	assert scores2.shape == (5, 20)
	assert offsets2.shape == (5, 20)
	assert overlaps2.shape == (5, 20)
	assert strands2.shape == (5, 20)

	p3, scores3, offsets3, overlaps3, strands3 = tomtom(pwms, pwms)

	assert p3.shape == (20, 20)
	assert scores3.shape == (20, 20)
	assert offsets3.shape == (20, 20)
	assert overlaps3.shape == (20, 20)
	assert strands3.shape == (20, 20)

	assert_array_almost_equal(p, p2[:2], 6)
	assert_array_almost_equal(scores, scores2[:2], 6)
	assert_array_almost_equal(offsets, offsets2[:2], 6)
	assert_array_almost_equal(overlaps, overlaps2[:2], 6)
	assert_array_almost_equal(strands, strands2[:2], 6)

	assert_array_almost_equal(p, p3[:2], 6)
	assert_array_almost_equal(scores, scores3[:2], 6)
	assert_array_almost_equal(offsets, offsets3[:2], 6)
	assert_array_almost_equal(overlaps, overlaps3[:2], 6)
	assert_array_almost_equal(strands, strands3[:2], 6)


def test_tomtom_self():
	pwms = generate_random_meme(n=1)
	p, scores, offsets, overlaps, strands = tomtom(pwms, pwms)

	assert p.shape == (1, 1)
	assert scores.shape == (1, 1)
	assert offsets.shape == (1, 1)
	assert overlaps.shape == (1, 1)
	assert strands.shape == (1, 1)

	assert_array_almost_equal(p, [[4.4409e-16]], 6)
	assert_array_almost_equal(scores, [[1346.]])
	assert_array_almost_equal(offsets, [[0.]])
	assert_array_almost_equal(overlaps, [[16.]])
	assert_array_almost_equal(strands, [[0.]])


def test_tomtom_selfp1():
	pwms = generate_random_meme(n=1)
	p, scores, offsets, overlaps, strands = tomtom(pwms, [p + 1 for p in pwms])

	assert p.shape == (1, 1)
	assert scores.shape == (1, 1)
	assert offsets.shape == (1, 1)
	assert overlaps.shape == (1, 1)
	assert strands.shape == (1, 1)

	assert_array_almost_equal(p, [[-1.7764e-15]], 6)
	assert_array_almost_equal(scores, [[1491.]])
	assert_array_almost_equal(offsets, [[0.]])
	assert_array_almost_equal(overlaps, [[16.]])
	assert_array_almost_equal(strands, [[0.]])


def test_tomtom_self_rc():
	pwms = generate_random_meme(n=1)
	p, scores, offsets, overlaps, strands = tomtom(pwms, 
		[p[::-1, ::-1] for p in pwms])

	assert p.shape == (1, 1)
	assert scores.shape == (1, 1)
	assert offsets.shape == (1, 1)
	assert overlaps.shape == (1, 1)
	assert strands.shape == (1, 1)

	assert_array_almost_equal(p, [[4.4409e-16]], 6)
	assert_array_almost_equal(scores, [[1346.]])
	assert_array_almost_equal(offsets, [[0.]])
	assert_array_almost_equal(overlaps, [[16.]])
	assert_array_almost_equal(strands, [[1.]])


def test_tomtom_meme():
	pwms = list(read_meme("tests/data/test.meme").values())
	p, scores, offsets, overlaps, strands = tomtom(pwms[:1], pwms)

	assert p.shape == (1, 12)
	assert scores.shape == (1, 12)
	assert offsets.shape == (1, 12)
	assert overlaps.shape == (1, 12)
	assert strands.shape == (1, 12)

	assert_array_almost_equal(p[0], [-1.687538e-14, 0.959270, 0.990233, 
		0.501984, 0.662968, 0.993437, 0.218161, 0.999998, 0.2650769, 0.53301, 
		0.872186, 0.71878], 4)
	assert_array_almost_equal(scores[0], [879.0, 557.0, 565.0, 617.0, 582.0, 
		573.0, 626.0, 515.0, 628.0, 607.0, 599.0, 587.0], 6)
	assert_array_almost_equal(offsets[0], [0.0, 0.0, 0.0, 7.0, -2.0, 2.0, 2.0, 
		2.0, 1.0, -2.0, 7.0, 0.0], 6)
	assert_array_almost_equal(overlaps[0], [10.0, 9.0, 10.0, 8.0, 8.0, 10.0, 
		8.0, 8.0, 10.0, 8.0, 10.0, 10.0], 6)
	assert_array_almost_equal(strands[0], [0., 1., 0., 1., 0., 1., 0., 0., 0., 
		1., 1., 0.])


def test_tomtom_reverse_complement():
	pwms = list(read_meme("tests/data/test.meme").values())
	p, scores, offsets, overlaps, strands = tomtom(pwms[:1], pwms, 
		reverse_complement=False)

	assert p.shape == (1, 12)
	assert scores.shape == (1, 12)
	assert offsets.shape == (1, 12)
	assert overlaps.shape == (1, 12)
	assert strands.shape == (1, 12)

	assert_array_almost_equal(p[0], [-0.0, 0.942776, 0.910807, 0.320109, 
		0.394581, 0.997576, 0.125284, 0.998832, 0.140479, 0.502672, 
		0.724997, 0.478826], 4)
	assert_array_almost_equal(scores[0], [877., 539., 563., 614., 584., 543., 
		624., 514., 628., 591., 592., 586.], 6)
	assert_array_almost_equal(offsets[0], [ 0., -1.,  0.,  6., -2.,  9.,  2.,  
		2.,  1.,  1.,  7.,  0.], 6)
	assert_array_almost_equal(overlaps[0], [10.,  9., 10.,  9.,  8., 10.,  8.,  
		8., 10., 10., 10., 10.], 6)
	assert_array_almost_equal(strands[0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


def test_tomtom_n_jobs():
	pwms = list(read_meme("tests/data/test.meme").values())
	p, scores, offsets, overlaps, strands = tomtom(pwms, pwms)
	p2, scores2, offsets2, overlaps2, strands2 = tomtom(pwms, pwms, n_jobs=1)

	assert_array_almost_equal(p, p2)
	assert_array_almost_equal(scores, scores2)
	assert_array_almost_equal(offsets, offsets2)
	assert_array_almost_equal(overlaps, overlaps2)
	assert_array_almost_equal(strands, strands2)


def test_tomtom_n_nearest():
	pwms = list(read_meme("tests/data/test.meme").values())
	p, scores, offsets, overlaps, strands, idxs = tomtom(pwms, pwms, 
		n_nearest=3)
	p2, scores2, offsets2, overlaps2, strands2 = tomtom(pwms, pwms)

	assert p.shape == (12, 3)
	assert scores.shape == (12, 3)
	assert offsets.shape == (12, 3)
	assert overlaps.shape == (12, 3)
	assert strands.shape == (12, 3)

	idxs = numpy.argsort(p2, axis=-1)[:, :3]
	for i, idx in enumerate(idxs):
		assert_array_almost_equal(p[i], p2[i, idx])
		assert_array_almost_equal(scores[i], scores2[i, idx])
		assert_array_almost_equal(offsets[i], offsets2[i, idx])
		assert_array_almost_equal(overlaps[i], overlaps2[i, idx])
		assert_array_almost_equal(strands[i], strands2[i, idx])


def test_tomtom_n_target_bins_small():
	pwms = list(read_meme("tests/data/test.meme").values())
	p, scores, offsets, overlaps, strands = tomtom(pwms[:1], pwms, 
		n_target_bins=10)

	assert p.shape == (1, 12)
	assert scores.shape == (1, 12)
	assert offsets.shape == (1, 12)
	assert overlaps.shape == (1, 12)
	assert strands.shape == (1, 12)

	assert_array_almost_equal(p[0], [-0.0, 0.921176, 0.994284, 0.525705, 
		0.653216, 0.985227, 0.214209, 0.999991, 0.304075, 0.455928, 0.848082, 
		0.788391], 4)
	assert_array_almost_equal(scores[0], [881., 565., 563., 617., 584., 580., 
		628., 521., 626., 614., 603., 583.], 6)
	assert_array_almost_equal(offsets[0], [ 0.,  0.,  0.,  6.,  0.,  2.,  2., 
		2.,  1., -2.,  7.,  0.], 6)
	assert_array_almost_equal(overlaps[0], [10.,  9., 10.,  9.,  8., 10.,  8.,  
		8., 10.,  8., 10., 10.], 6)
	assert_array_almost_equal(strands[0], [0., 1., 0., 0., 1., 1., 0., 0., 0.,
		1., 1., 0.])


###


def test_binned_median_single():
	X = numpy.array([5.0], dtype='float64')
	bins = numpy.zeros((3, 2), dtype='float64')
	counts = numpy.ones(1, dtype='int64')
	assert _binned_median(X, bins, 5.0, 6.0, counts) == 5.0


def test_binned_median_identical():
	X = numpy.array([3.0, 3.0, 3.0, 3.0], dtype='float64')
	bins = numpy.zeros((5, 2), dtype='float64')
	counts = numpy.ones(4, dtype='int64')
	assert _binned_median(X, bins, 3.0, 4.0, counts) == 3.0


def test_binned_median_dominant_bin():
	X = numpy.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype='float64')
	bins = numpy.zeros((5, 2), dtype='float64')
	counts = numpy.ones(5, dtype='int64')
	counts[0] = 100
	assert _binned_median(X, bins, 0.0, 4.0, counts) == 0.0


###


def test_tomtom_n_target_bins_none():
	pwms = list(read_meme("tests/data/test.meme").values())
	p, scores, offsets, overlaps, strands = tomtom(pwms[:1], pwms,
		n_target_bins=None)

	assert p.shape == (1, 12)
	assert scores.shape == (1, 12)
	assert offsets.shape == (1, 12)
	assert overlaps.shape == (1, 12)
	assert strands.shape == (1, 12)

	assert_array_almost_equal(p[0], [-2.087219e-14, 9.594323e-01,
		9.903021e-01, 5.028105e-01, 6.634917e-01, 9.923712e-01, 2.186044e-01,
		9.999984e-01, 2.656418e-01, 5.193330e-01, 8.727407e-01, 7.193744e-01],
		4)
	assert_array_almost_equal(scores[0], [879., 557., 565., 617., 582., 574.,
		626., 515., 628., 608., 599., 587.], 6)
	assert_array_almost_equal(offsets[0], [0., 0., 0., 7., -2., 2., 2., 2., 1.,
		-2., 7., 0.], 6)
	assert_array_almost_equal(overlaps[0], [10., 9., 10., 8., 8., 10., 8., 8.,
		10., 8., 10., 10.], 6)
	assert_array_almost_equal(strands[0], [0., 1., 0., 1., 0., 1., 0., 0., 0.,
		1., 1., 0.])


def test_tomtom_n_target_bins_none_vs_default():
	pwms = list(read_meme("tests/data/test.meme").values())
	p0, s0, o0, v0, t0 = tomtom(pwms[:1], pwms, n_target_bins=None)
	p1, s1, o1, v1, t1 = tomtom(pwms[:1], pwms)

	# Hashing is approximate, so the no-hashing path is allowed to differ
	# from the default but should be in the same ballpark.
	assert numpy.abs(p0 - p1).max() < 0.05
	assert numpy.abs(s0 - s1).max() <= 10


def test_tomtom_n_nearest_one():
	pwms = list(read_meme("tests/data/test.meme").values())
	p, scores, offsets, overlaps, strands, idxs = tomtom(pwms, pwms,
		n_nearest=1)
	p2, scores2, offsets2, overlaps2, strands2 = tomtom(pwms, pwms)

	assert p.shape == (12, 1)
	assert scores.shape == (12, 1)
	assert offsets.shape == (12, 1)
	assert overlaps.shape == (12, 1)
	assert strands.shape == (12, 1)
	assert idxs.shape == (12, 1)

	for i in range(len(pwms)):
		idx = idxs[i].astype(int)
		assert_array_almost_equal(p[i], p2[i, idx])
		assert_array_almost_equal(scores[i], scores2[i, idx])
		assert_array_almost_equal(offsets[i], offsets2[i, idx])
		assert_array_almost_equal(overlaps[i], overlaps2[i, idx])
		assert_array_almost_equal(strands[i], strands2[i, idx])


def test_tomtom_n_nearest_all():
	pwms = list(read_meme("tests/data/test.meme").values())
	K = len(pwms)
	p, scores, offsets, overlaps, strands, idxs = tomtom(pwms, pwms,
		n_nearest=K)
	p2, scores2, offsets2, overlaps2, strands2 = tomtom(pwms, pwms)

	assert p.shape == (12, K)
	assert idxs.shape == (12, K)

	for i in range(len(pwms)):
		idx = idxs[i].astype(int)
		assert_array_almost_equal(p[i], p2[i, idx])
		assert_array_almost_equal(scores[i], scores2[i, idx])
		assert_array_almost_equal(offsets[i], offsets2[i, idx])
		assert_array_almost_equal(overlaps[i], overlaps2[i, idx])
		assert_array_almost_equal(strands[i], strands2[i, idx])


def test_tomtom_n_nearest_sorted():
	pwms = list(read_meme("tests/data/test.meme").values())
	for K in (1, 3, len(pwms)):
		p, scores, offsets, overlaps, strands, idxs = tomtom(pwms, pwms,
			n_nearest=K)

		# p-values must be returned in non-decreasing order per query
		for i in range(len(pwms)):
			assert numpy.all(numpy.diff(p[i]) >= -1e-9)


def test_tomtom_different_lengths_and_counts():
	pwms = generate_random_meme(n=20)
	p, scores, offsets, overlaps, strands = tomtom(pwms[:3], pwms[5:9])

	assert p.shape == (3, 4)
	assert scores.shape == (3, 4)
	assert offsets.shape == (3, 4)
	assert overlaps.shape == (3, 4)
	assert strands.shape == (3, 4)

	# Selecting a subset of queries (with the same target set) must give
	# identical results to running the full query set.
	pf, sf, of, vf, tf = tomtom(pwms, pwms[5:9])

	assert_array_almost_equal(p, pf[:3], 6)
	assert_array_almost_equal(scores, sf[:3], 6)
	assert_array_almost_equal(offsets, of[:3], 6)
	assert_array_almost_equal(overlaps, vf[:3], 6)
	assert_array_almost_equal(strands, tf[:3], 6)


def test_tomtom_single_column():
	state = numpy.random.RandomState(0)
	single = []
	for i in range(4):
		pwm = state.rand(4, 1)
		single.append(pwm / pwm.sum(axis=0, keepdims=True))

	p, scores, offsets, overlaps, strands = tomtom(single, single)

	assert p.shape == (4, 4)

	assert_array_almost_equal(p, [
		[0.234375, 0.609375, 0.75, 0.984375],
		[0.609375, 0.234375, 0.4375, 0.984375],
		[0.609375, 0.4375, 0.234375, 0.984375],
		[0.4375, 0.859375, 0.609375, 0.234375]], 4)
	assert_array_almost_equal(scores, [
		[99., 85., 84., 66.],
		[86., 99., 95., 64.],
		[86., 95., 99., 66.],
		[72., 68., 69., 99.]], 6)
	assert_array_almost_equal(overlaps, numpy.ones((4, 4)), 6)
	assert_array_almost_equal(offsets, numpy.zeros((4, 4)), 6)


def test_tomtom_mixed_short_long():
	state = numpy.random.RandomState(2)

	def mk(length):
		pwm = state.rand(4, length)
		return pwm / pwm.sum(axis=0, keepdims=True)

	mixed = [mk(1), mk(2), mk(3), mk(25), mk(30)]
	p, scores, offsets, overlaps, strands = tomtom(mixed, mixed)

	assert p.shape == (5, 5)
	assert scores.shape == (5, 5)
	assert offsets.shape == (5, 5)
	assert overlaps.shape == (5, 5)
	assert strands.shape == (5, 5)

	# Each motif matches itself best (lowest p-value on the diagonal).
	assert_array_almost_equal(numpy.diag(p), [1.632626e-02, 1.343680e-04,
		1.101413e-06, -7.149836e-14, -1.039169e-13], 6)

	# The overlap of a query against itself is the full motif length.
	assert_array_almost_equal(numpy.diag(overlaps), [1., 2., 3., 25., 30.], 6)
	assert_array_almost_equal(numpy.diag(offsets), [0., 0., 0., 0., 0.], 6)

	# Overlap with a target can never exceed the shorter motif's length.
	assert_array_almost_equal(overlaps, [
		[1., 1., 1., 1., 1.],
		[1., 2., 1., 2., 2.],
		[1., 2., 3., 3., 3.],
		[1., 2., 3., 25., 25.],
		[1., 2., 3., 25., 30.]], 6)


def test_tomtom_reverse_complement_merge():
	# An end-to-end check of the reverse-complement merge that complements the
	# helper-level `test_merge_rc_results`. Note that the per-strand p-values
	# and scores from an RC run cannot be reconstructed exactly by two separate
	# `reverse_complement=False` runs, because the score binning, medians, and
	# background distributions are estimated over the *combined* target set
	# (forward + reverse) in a single RC run. We therefore check the properties
	# that hold exactly within the RC run plus golden values.
	pwms = list(read_meme("tests/data/test.meme").values())

	p_rc, s_rc, o_rc, v_rc, t_rc = tomtom(pwms[:1], pwms,
		reverse_complement=True)

	# Run both strands separately with reverse_complement=False, which gives
	# the strand-selection direction (which strand scores higher).
	p_f, s_f, o_f, v_f, t_f = tomtom(pwms[:1], pwms, reverse_complement=False)
	p_r, s_r, o_r, v_r, t_r = tomtom(pwms[:1],
		[t[::-1, ::-1] for t in pwms], reverse_complement=False)

	# Strands are binary.
	assert set(numpy.unique(t_rc)).issubset({0.0, 1.0})

	# The RC run selects the reverse strand exactly where the reverse strand
	# scores strictly higher than the forward strand.
	rev_wins = (s_r > s_f).astype('float64')
	assert_array_almost_equal(t_rc, rev_wins)

	# Golden values captured from the RC run.
	assert_array_almost_equal(p_rc[0], [-1.687538e-14, 0.959270, 0.990233,
		0.501984, 0.662968, 0.993437, 0.218161, 0.999998, 0.265077, 0.533010,
		0.872186, 0.718780], 4)
	assert_array_almost_equal(s_rc[0], [879., 557., 565., 617., 582., 573.,
		626., 515., 628., 607., 599., 587.], 6)
	assert_array_almost_equal(t_rc[0], [0., 1., 0., 1., 0., 1., 0., 0., 0., 1.,
		1., 0.])


def test_tomtom_p_values_non_negative():
	# Regression test for the small, negative p-values that users reported on
	# very good matches (#7).
	#
	# The p-value of a hit is read out of the background survival function
	# `B` built in `_p_value_backgrounds`. When it was formed as
	# `1 - cumsum(pdf)`, round-off in the cumsum over thousands of bins made
	# the survival value of the very best matches, in the extreme right tail
	# where the CDF is ~1, a tiny negative number (~ -1e-14). It is now the
	# sum of the pdf above each score, which is never negative, and clamped
	# to [0, 1].
	#
	# A self-comparison of `test.meme` exercises this: every motif's best hit
	# is itself, sitting in that tail. The p-values must be non-negative on
	# both strands; `test_tomtom_self_matches_positive` checks that they are
	# also not rounded to 0.
	pwms = list(read_meme("tests/data/test.meme").values())

	for rc in (True, False):
		p, scores, offsets, overlaps, strands = tomtom(pwms, pwms,
			reverse_complement=rc)
		assert numpy.all(p >= 0), p.min()


def test_tomtom_n_jobs_subsets():
	pwms = list(read_meme("tests/data/test.meme").values())
	out1 = tomtom(pwms[:5], pwms, n_jobs=1)
	out2 = tomtom(pwms[:5], pwms, n_jobs=-1)

	for a, b in zip(out1, out2):
		assert_array_almost_equal(a, b)


def test_tomtom_n_jobs_n_nearest():
	pwms = list(read_meme("tests/data/test.meme").values())
	out1 = tomtom(pwms, pwms, n_nearest=4, n_jobs=1)
	out2 = tomtom(pwms, pwms, n_nearest=4, n_jobs=-1)

	for a, b in zip(out1, out2):
		assert_array_almost_equal(a, b)


##
# Helpers for the property tests below.


def _random_pwms(lengths, random_state=0):
	"""Random PWMs with the given lengths, built like generate_random_meme."""

	state = numpy.random.RandomState(random_state)

	pwms = []
	for length in lengths:
		pwm = state.choice(17, p=[0.2] + [0.05]*16, size=(length, 4))
		pwm[pwm.sum(axis=1) == 0] = [1, 0, 0, 0]
		pwm = pwm / pwm.sum(axis=1, keepdims=True)
		pwms.append(pwm.T)

	return pwms


def _lengths(regime, n, is_query, random_state):
	state = numpy.random.RandomState(random_state)

	if regime == 'short':
		return state.randint(1, 6, size=n)
	elif regime == 'equal':
		return numpy.full(n, 8)
	elif regime == 'mixed':
		return state.randint(1, 33, size=n)
	elif regime == 'qlong':
		return state.randint(20, 31, size=n) if is_query else \
			state.randint(3, 9, size=n)
	elif regime == 'tlong':
		return state.randint(3, 9, size=n) if is_query else \
			state.randint(20, 31, size=n)


def _assert_identical(out0, out1):
	assert len(out0) == len(out1)
	for a, b in zip(out0, out1):
		assert a.shape == b.shape
		assert_array_equal(a, b)


def _assert_valid(out, Qs, Ts, n_nearest=None, reverse_complement=True):
	"""Structural checks that any correct TOMTOM output must satisfy."""

	n_q, n_t = len(Qs), len(Ts)
	n_out = n_t if n_nearest is None else n_nearest

	assert isinstance(out, numpy.ndarray)
	assert len(out) == (5 if n_nearest is None else 6)
	for x in out:
		assert x.dtype == numpy.float64
		assert x.shape == (n_q, n_out)
		assert numpy.all(numpy.isfinite(x))

	p, scores, offsets, overlaps, strands = out[:5]

	if n_nearest is None:
		t_idxs = numpy.tile(numpy.arange(n_t), (n_q, 1))
	else:
		t_idxs = out[5].astype(int)
		assert_array_equal(out[5], t_idxs)
		assert t_idxs.min() >= 0 and t_idxs.max() < n_t
		for row in t_idxs:
			assert len(numpy.unique(row)) == len(row)

	assert p.min() >= 0 and p.max() <= 1
	assert set(numpy.unique(strands)).issubset({0.0, 1.0})
	if not reverse_complement:
		assert numpy.all(strands == 0)

	for x in (scores, offsets, overlaps):
		assert_array_equal(x, numpy.round(x))

	nq = numpy.array([Q.shape[-1] for Q in Qs])[:, None]
	nt = numpy.array([T.shape[-1] for T in Ts])[t_idxs]

	assert numpy.all(overlaps >= 1)
	assert numpy.all(overlaps <= numpy.minimum(nq, nt))
	assert numpy.all(offsets >= -(nq - 1))
	assert numpy.all(offsets <= nt - 1)

	k = offsets + nq - 1
	assert_array_equal(overlaps, numpy.minimum(k + 1, nq) -
		numpy.maximum(0, k - nt + 1))


##


@pytest.mark.parametrize("n_q", [1, 2, 7, 33])
@pytest.mark.parametrize("n_t", [1, 2, 9, 40])
@pytest.mark.parametrize("regime", ['short', 'equal', 'mixed', 'qlong', 
	'tlong'])
def test_tomtom_shape_grid(n_q, n_t, regime):
	Qs = _random_pwms(_lengths(regime, n_q, True, n_q), random_state=n_q)
	Ts = _random_pwms(_lengths(regime, n_t, False, 100+n_t), 
		random_state=100+n_t)

	out = tomtom(Qs, Ts, n_jobs=_N_JOBS)
	_assert_valid(out, Qs, Ts)


@pytest.mark.parametrize("n_nearest", [1, 2, 5, 9])
@pytest.mark.parametrize("regime", ['short', 'equal', 'mixed', 'qlong', 
	'tlong'])
def test_tomtom_shape_grid_n_nearest(n_nearest, regime):
	Qs = _random_pwms(_lengths(regime, 7, True, 0), random_state=0)
	Ts = _random_pwms(_lengths(regime, 9, False, 1), random_state=1)

	out = tomtom(Qs, Ts, n_nearest=n_nearest, n_jobs=_N_JOBS)
	_assert_valid(out, Qs, Ts, n_nearest=n_nearest)

	full = tomtom(Qs, Ts, n_jobs=_N_JOBS)
	idxs = out[5].astype(int)
	for i in range(len(Qs)):
		assert numpy.all(numpy.diff(out[0][i]) >= 0)
		for x, y in zip(out[:5], full):
			assert_array_equal(x[i], y[i, idxs[i]])

		# The kept targets are the n_nearest smallest p-values.
		assert_array_equal(out[0][i], numpy.sort(full[0][i])[:n_nearest])


@pytest.mark.parametrize("regime", ['short', 'equal', 'mixed', 'qlong', 
	'tlong'])
def test_tomtom_shape_grid_no_rc(regime):
	Qs = _random_pwms(_lengths(regime, 7, True, 0), random_state=0)
	Ts = _random_pwms(_lengths(regime, 9, False, 1), random_state=1)

	out = tomtom(Qs, Ts, reverse_complement=False, n_jobs=_N_JOBS)
	_assert_valid(out, Qs, Ts, reverse_complement=False)


def test_tomtom_n_nearest_larger_than_targets():
	Qs = _random_pwms([6, 8, 10], random_state=0)
	Ts = _random_pwms([5, 7, 9, 11, 4], random_state=1)

	try:
		out = tomtom(Qs, Ts, n_nearest=8)
	except ValueError:
		return

	_assert_valid(out, Qs, Ts, n_nearest=out[0].shape[1])
	assert out[0].shape[1] <= len(Ts)


##


@pytest.fixture
def mixed_pwms():
	lengths = numpy.random.RandomState(3).randint(1, 31, size=20)
	return _random_pwms(lengths, random_state=3)


def test_tomtom_query_permutation(mixed_pwms):
	out = tomtom(mixed_pwms, mixed_pwms)

	for seed in range(3):
		perm = numpy.random.RandomState(seed).permutation(len(mixed_pwms))
		out_p = tomtom([mixed_pwms[i] for i in perm], mixed_pwms)
		_assert_identical([x[perm] for x in out], out_p)


def test_tomtom_query_duplicates(mixed_pwms):
	out = tomtom([mixed_pwms[3]]*4 + [mixed_pwms[5]]*3, mixed_pwms)

	for x in out:
		for k in range(1, 4):
			assert_array_equal(x[0], x[k])
		assert_array_equal(x[4], x[5])
		assert_array_equal(x[4], x[6])


def test_tomtom_query_subsets(mixed_pwms):
	out = tomtom(mixed_pwms, mixed_pwms)
	state = numpy.random.RandomState(0)

	for _ in range(10):
		n = state.randint(1, len(mixed_pwms))
		idxs = numpy.sort(state.choice(len(mixed_pwms), size=n, replace=False))
		out_s = tomtom([mixed_pwms[i] for i in idxs], mixed_pwms)
		_assert_identical([x[idxs] for x in out], out_s)


@pytest.mark.parametrize("reverse_complement", [True, False])
def test_tomtom_target_permutation(mixed_pwms, reverse_complement):
	# Exact only without target hashing: with hashing, the representative
	# column kept for each hash bin is the first one seen, so the target
	# order changes the (approximate) columns that are scored.
	out = tomtom(mixed_pwms, mixed_pwms, n_target_bins=None, 
		reverse_complement=reverse_complement)

	for seed in range(3):
		perm = numpy.random.RandomState(seed).permutation(len(mixed_pwms))
		out_p = tomtom(mixed_pwms, [mixed_pwms[i] for i in perm], 
			n_target_bins=None, reverse_complement=reverse_complement)
		_assert_identical([x[:, perm] for x in out], out_p)


@pytest.mark.parametrize("n_jobs", [1, 2, 3, 5, -1])
@pytest.mark.parametrize("n_nearest", [None, 3])
def test_tomtom_n_jobs_identical(mixed_pwms, n_jobs, n_nearest):
	_require_threads(n_jobs)
	out = tomtom(mixed_pwms, mixed_pwms, n_nearest=n_nearest, n_jobs=1)
	out_j = tomtom(mixed_pwms, mixed_pwms, n_nearest=n_nearest, n_jobs=n_jobs)
	_assert_identical(out, out_j)


@pytest.mark.parametrize("n_jobs", [1, 2, 3, -1])
def test_tomtom_n_jobs_restores_threads(mixed_pwms, n_jobs):
	_require_threads(n_jobs)
	before = numba.get_num_threads()
	tomtom(mixed_pwms[:3], mixed_pwms, n_jobs=n_jobs)
	assert numba.get_num_threads() == before


def test_tomtom_inputs_not_mutated(mixed_pwms):
	Qs = [x.copy() for x in mixed_pwms[:5]]
	Ts = [x.copy() for x in mixed_pwms]
	Qs_list, Ts_list = list(Qs), list(Ts)

	tomtom(Qs, Ts)
	tomtom(Qs, Ts, reverse_complement=False, n_target_bins=None)

	assert len(Qs) == 5 and len(Ts) == len(mixed_pwms)
	for a, b in zip(Qs, mixed_pwms[:5]):
		assert_array_equal(a, b)
	for a, b in zip(Ts, mixed_pwms):
		assert_array_equal(a, b)
	assert all(a is b for a, b in zip(Qs, Qs_list))
	assert all(a is b for a, b in zip(Ts, Ts_list))


def test_tomtom_non_contiguous(mixed_pwms):
	out = tomtom(mixed_pwms, mixed_pwms)

	Qs = [numpy.asfortranarray(x) for x in mixed_pwms]
	Ts = [numpy.repeat(x, 2, axis=1)[:, ::2] for x in mixed_pwms]
	assert not Ts[0].flags['C_CONTIGUOUS'] or Ts[0].shape[1] == 1
	_assert_identical(out, tomtom(Qs, Ts))

	# Reversed views of reversed arrays are non-contiguous views of the
	# original values.
	rr = [numpy.ascontiguousarray(x[::-1, ::-1])[::-1, ::-1] 
		for x in mixed_pwms]
	_assert_identical(out, tomtom(rr, rr))


def test_tomtom_float32(mixed_pwms):
	# float32 inputs are not converted to float64, so the distances are 
	# computed at lower precision and a few scores move by one bin. The
	# results must stay valid and agree closely with float64.
	out = tomtom(mixed_pwms, mixed_pwms)

	Xs = [x.astype('float32') for x in mixed_pwms]
	out32 = tomtom(Xs, Xs)

	_assert_valid(out32, Xs, Xs)
	assert numpy.abs(out[1] - out32[1]).max() <= 1
	assert (out[1] != out32[1]).mean() < 0.05
	assert_array_equal(out[0].argmin(axis=1), out32[0].argmin(axis=1))


def test_tomtom_query_with_zero_column(mixed_pwms):
	# A single all-zero column is allowed; only all-zero inputs are rejected.
	Q = mixed_pwms[0].copy()
	Q[:, 0] = 0
	Q = numpy.concatenate([Q, mixed_pwms[1]], axis=-1)

	out = tomtom([Q], mixed_pwms)
	_assert_valid(out, [Q], mixed_pwms)


@pytest.mark.parametrize("n_q", [1, 3])
def test_tomtom_all_zero_raises(mixed_pwms, n_q):
	zeros = numpy.zeros((4, 7))

	assert_raises(ValueError, tomtom, [zeros]*n_q, mixed_pwms)
	assert_raises(ValueError, tomtom, mixed_pwms, [zeros]*n_q)
	assert_raises(ValueError, tomtom, mixed_pwms, [zeros]*n_q, 
		reverse_complement=False)
	assert_raises(ValueError, tomtom, [zeros], [zeros], n_target_bins=None)


##


def test_tomtom_reverse_complement_targets_untied():
	# Reversing every target swaps the strands of the best alignment. Pairs
	# whose forward and reverse best scores tie are resolved to strand 1 in
	# both runs (see _merge_rc_results), so they are excluded.
	pwms = generate_random_meme(n=20)
	rc = [p[::-1, ::-1] for p in pwms]

	out0 = tomtom(pwms, pwms)
	out1 = tomtom(pwms, rc)

	assert_array_equal(out0[0], out1[0])
	assert_array_equal(out0[1], out1[1])

	s_f = tomtom(pwms, pwms, reverse_complement=False)[1]
	s_r = tomtom(pwms, rc, reverse_complement=False)[1]

	untied = (out0[4] + out1[4]) == 1
	assert untied.mean() > 0.9

	for x0, x1 in zip(out0[2:4], out1[2:4]):
		assert_array_equal(x0[untied], x1[untied])

	# Where strands are tied, both runs report strand 1.
	assert_array_equal(out0[4][~untied], 1)
	assert_array_equal(out1[4][~untied], 1)


def test_tomtom_self_is_best(mixed_pwms):
	informative = [x for x in mixed_pwms if x.shape[-1] >= 6]
	p, scores, offsets, overlaps, strands = tomtom(informative, informative)

	n = len(informative)
	assert_array_equal(p.argmin(axis=1), numpy.arange(n))
	assert_array_equal(numpy.diag(offsets), numpy.zeros(n))
	assert_array_equal(numpy.diag(overlaps), 
		[x.shape[-1] for x in informative])
	assert_array_equal(numpy.diag(strands), numpy.zeros(n))
	assert numpy.all(numpy.diag(p) < 1e-4)

	# A PWM's score against itself is the largest in its row.
	assert_array_equal(scores.argmax(axis=1), numpy.arange(n))


@pytest.mark.parametrize("j", range(12))
def test_tomtom_substring_queries(j):
	# Every substring of a target, and its reverse complement, is placed at
	# the position it was cut from with the correct strand.
	pwms = list(read_meme("tests/data/test.meme").values())
	t = pwms[j]
	L = t.shape[-1]

	Qs, spans = [], []
	for a in range(0, L-3):
		for b in range(a+4, L+1):
			Qs.append(t[:, a:b])
			Qs.append(t[:, a:b][::-1, ::-1])
			spans.append((a, b))

	p, scores, offsets, overlaps, strands = tomtom(Qs, pwms)

	for i, (a, b) in enumerate(spans):
		f, r = 2*i, 2*i + 1

		assert p[f].argmin() == j
		assert offsets[f, j] == a
		assert overlaps[f, j] == b - a
		assert strands[f, j] == 0

		assert p[r].argmin() == j
		assert offsets[r, j] == L - b
		assert overlaps[r, j] == b - a
		assert strands[r, j] == 1


def test_tomtom_substring_queries_no_rc():
	pwms = list(read_meme("tests/data/test.meme").values())

	for j, t in enumerate(pwms):
		L = t.shape[-1]
		Qs = [t[:, a:a+5] for a in range(L-4)]
		p, scores, offsets, overlaps, strands = tomtom(Qs, pwms, 
			reverse_complement=False)

		assert_array_equal(p.argmin(axis=1), j)
		assert_array_equal(offsets[:, j], numpy.arange(L-4))
		assert_array_equal(overlaps[:, j], 5)
		assert_array_equal(strands, 0)


def test_tomtom_rc_matches_query_rc(mixed_pwms):
	# Reverse complementing the targets under reverse_complement=False scores
	# the same set of alignments as reverse complementing the queries, so the
	# best scores match exactly and the p-values to round-off. Offsets and
	# overlaps are not compared because equal-scoring alignments are broken
	# in scan order, which differs between the two orientations.
	Qs = mixed_pwms[:6]
	Ts = mixed_pwms

	out = tomtom(Qs, [T[::-1, ::-1] for T in Ts], reverse_complement=False,
		n_target_bins=None)
	out_q = tomtom([Q[::-1, ::-1] for Q in Qs], Ts, reverse_complement=False,
		n_target_bins=None)

	_assert_valid(out, Qs, Ts, reverse_complement=False)
	assert_array_equal(out[1], out_q[1])
	assert_array_almost_equal(out[0], out_q[0], 12)


##


# 255 and 400 push the histogram offset past the default n_cache and, at 400,
# past the range of the old int8 `_gamma_int`.
@pytest.mark.parametrize("n_score_bins", [10, 50, 100, 127, 150, 200, 255, 
	400])
def test_tomtom_n_score_bins(n_score_bins):
	pwms = list(read_meme("tests/data/test.meme").values())
	out = tomtom(pwms, pwms, n_score_bins=n_score_bins)

	_assert_valid(out, pwms, pwms)
	p, scores = out[:2]
	assert_array_equal(p.argmin(axis=1), numpy.arange(12))
	assert_array_equal(scores.argmax(axis=1), numpy.arange(12))

	# Scores are sums of per-column bins, so they scale with the number of
	# bins.
	ref = tomtom(pwms, pwms)[1]
	ratio = scores / ref
	assert numpy.abs(numpy.median(ratio) - n_score_bins / 100.) < 0.1


@pytest.mark.parametrize("n_median_bins", [10, 100, 1000, 5000])
def test_tomtom_n_median_bins(n_median_bins):
	pwms = list(read_meme("tests/data/test.meme").values())
	out = tomtom(pwms, pwms, n_median_bins=n_median_bins)

	_assert_valid(out, pwms, pwms)
	assert_array_equal(out[0].argmin(axis=1), numpy.arange(12))

	ref = tomtom(pwms, pwms, n_median_bins=5000)
	if n_median_bins >= 1000:
		assert numpy.abs(out[0] - ref[0]).max() < 0.05


@pytest.mark.parametrize("n_target_bins", [None, 3, 5, 10, 100, 1000])
def test_tomtom_n_target_bins(n_target_bins, mixed_pwms):
	pwms = list(read_meme("tests/data/test.meme").values())
	out = tomtom(pwms, pwms, n_target_bins=n_target_bins)

	_assert_valid(out, pwms, pwms)
	if n_target_bins is None or n_target_bins >= 10:
		assert_array_equal(out[0].argmin(axis=1), numpy.arange(12))

	out = tomtom(mixed_pwms, mixed_pwms, n_target_bins=n_target_bins)
	_assert_valid(out, mixed_pwms, mixed_pwms)


def test_tomtom_n_target_bins_2(mixed_pwms):
	# Regression: with coarse hashing there are fewer unique target columns
	# than the longest target, which overran `t_sums` in `_p_values`.
	out = tomtom(mixed_pwms, mixed_pwms, n_target_bins=2)
	_assert_valid(out, mixed_pwms, mixed_pwms)


def test_p_values_fewer_unique_columns_than_target_length():
	# Regression: `t_sums` was sized by the number of unique target columns,
	# so a target longer than that overran it. The overrun corrupts the heap
	# without changing the output, so it is only visible to a bounds-checked
	# copy of the kernel.
	_p_values_checked = numba.njit(boundscheck=True)(_p_values.py_func)

	state = numpy.random.RandomState(0)
	nq, offset, n_bins = 3, 2, 10
	T_lens = numpy.array([12, 4], dtype='int64')
	gamma_int = state.randint(-offset, n_bins - offset + 1, 
		size=(2, nq)).astype('int16')
	rr_inv = state.randint(0, 2, size=T_lens.sum()).astype('uint64')
	B_cdfs = numpy.linspace(1, 0, nq*(n_bins+offset))[None].repeat(13, axis=0)
	results = numpy.zeros((2, 5))

	_p_values_checked(gamma_int, B_cdfs, rr_inv, T_lens, -1, nq, offset, 
		results)

	# The best score is the best sum along any alignment of the query.
	start = 0
	for i, nt in enumerate(T_lens):
		cols = gamma_int[rr_inv[start:start+nt].astype(int)]
		sums = numpy.full(nt+nq-1, nq*offset)
		for k in range(nt):
			sums[k:k+nq] += cols[k]

		assert results[i, 1] == sums.max()
		start += nt


def test_tomtom_n_target_bins_many_equals_none():
	# With a very fine hash every distinct column gets its own bin, so the
	# hashed and unhashed paths score the same columns.
	pwms = list(read_meme("tests/data/test.meme").values())
	out0 = tomtom(pwms, pwms, n_target_bins=None)
	out1 = tomtom(pwms, pwms, n_target_bins=100000)

	for a, b in zip(out0[1:], out1[1:]):
		assert_array_equal(a, b)
	assert_array_almost_equal(out0[0], out1[0], 10)


@pytest.mark.parametrize("n_cache", [100, 250, 500])
def test_tomtom_n_cache(n_cache):
	# n_cache only sizes the scratchpad, so it must not change the results.
	pwms = list(read_meme("tests/data/test.meme").values())
	_assert_identical(tomtom(pwms, pwms), tomtom(pwms, pwms, n_cache=n_cache))


@pytest.mark.parametrize("n_cache", [0, 5, 20])
def test_tomtom_n_cache_too_small(mixed_pwms, n_cache):
	# Regression: an offset above n_cache used to overrun the workspace. Such
	# queries now get their own workspace, so the results are unchanged.
	_assert_identical(tomtom(mixed_pwms[:5], mixed_pwms), 
		tomtom(mixed_pwms[:5], mixed_pwms, n_cache=n_cache))


##


@pytest.mark.parametrize("seed", range(5))
def test_binned_median_random_counts(seed):
	state = numpy.random.RandomState(seed)
	n, n_bins = state.randint(5, 2000), 1000

	X = state.randn(n)
	counts = state.randint(1, 10, size=n)
	bins = numpy.zeros((n_bins, 2), dtype='float64')

	m = _binned_median(X, bins, X.min(), X.max(), counts)
	expected = numpy.median(numpy.repeat(X, counts))

	width = (X.max() - X.min()) / (n_bins - 1)
	assert abs(m - expected) <= width
	assert X.min() <= m <= X.max()


@pytest.mark.parametrize("n_bins", [2, 10, 100, 10000])
def test_binned_median_bin_count(n_bins):
	X = numpy.random.RandomState(0).uniform(-5, 0, size=501)
	counts = numpy.ones(len(X), dtype='int64')
	bins = numpy.zeros((n_bins, 2), dtype='float64')

	m = _binned_median(X, bins, X.min(), X.max(), counts)
	width = (X.max() - X.min()) / (n_bins - 1)
	assert abs(m - numpy.median(X)) <= width


def test_binned_median_reuses_bins():
	# Bins are cleared on every call, so dirty scratch space is harmless.
	X = numpy.array([0, 4, 2, 1, 3], dtype='float64')
	counts = numpy.ones(5, dtype='int64')
	bins = numpy.full((5, 2), 1000.0)

	assert _binned_median(X, bins, 0, 4, counts) == 2
	assert _binned_median(X, bins, 0, 4, counts) == 2


@pytest.mark.parametrize("seed", range(4))
def test_pairwise_max_brute_force(seed):
	state = numpy.random.RandomState(seed)
	n = state.randint(1, 60)

	x = state.dirichlet(numpy.ones(n) * 0.5)
	y = state.dirichlet(numpy.ones(n) * 0.5)

	expected = numpy.zeros(n)
	for i in range(n):
		for j in range(n):
			expected[max(i, j)] += x[i] * y[j]

	z = numpy.empty(n)
	_pairwise_max(x, y, numpy.cumsum(y), z, n)
	assert_array_almost_equal(z, expected, 12)
	assert_array_almost_equal([z.sum()], [1.0], 12)


def test_pairwise_max_in_place():
	# `_p_value_backgrounds` passes the same array as `x` and `z`.
	state = numpy.random.RandomState(0)
	x = state.dirichlet(numpy.ones(20))
	y = state.dirichlet(numpy.ones(20))

	expected = numpy.empty(20)
	_pairwise_max(x.copy(), y, numpy.cumsum(y), expected, 20)

	_pairwise_max(x, y, numpy.cumsum(y), x, 20)
	assert_array_almost_equal(x, expected, 12)


def test_pairwise_max_partial_n():
	x = numpy.random.RandomState(0).dirichlet(numpy.ones(10))
	y = numpy.random.RandomState(1).dirichlet(numpy.ones(10))

	z = numpy.full(10, -7.0)
	_pairwise_max(x, y, numpy.cumsum(y), z, 6)

	assert_array_equal(z[6:], -7.0)
	assert_array_almost_equal(z[:6], (x * numpy.cumsum(y) + y * numpy.cumsum(x)
		- x * y)[:6], 12)


def test_merge_rc_results_ties():
	# Equal scores resolve to the reverse strand.
	results = numpy.array([
		[0.1, 5, 1, 3, 0],
		[0.2, 7, 2, 4, 0],
		[0.3, 9, 3, 5, 0],
		[0.4, 5, -1, 2, 0],
		[0.5, 6, -2, 1, 0],
		[0.1, 9, -3, 6, 0]
	], dtype='float64')

	_merge_rc_results(results)

	assert_array_almost_equal(results[:3, 0], [1 - 0.9**2, 1 - 0.8**2, 
		1 - 0.9**2], 12)
	assert_array_equal(results[:3, 1], [5, 7, 9])
	assert_array_equal(results[:3, 2], [-1, 2, -3])
	assert_array_equal(results[:3, 3], [2, 4, 6])
	assert_array_equal(results[:3, 4], [1, 0, 1])

	# The reverse-strand half is left untouched.
	assert_array_equal(results[3:, 1], [5, 6, 9])


@pytest.mark.parametrize("nq", [1, 2, 3, 6])
@pytest.mark.parametrize("t_max", [1, 3, 8])
def test_p_value_backgrounds_survival(nq, t_max):
	# B[t] is a survival function over integer scores for targets of length
	# t: it starts near 1, never increases, and stays within [0, 1].
	n_bins, n_cache, offset = 20, 20, 5
	n_len = nq * n_bins + nq * n_cache

	state = numpy.random.RandomState(nq * 10 + t_max)
	f = numpy.zeros((nq, n_bins + 1))
	f[:, 1:] = state.dirichlet(numpy.ones(n_bins), size=nq)

	A = numpy.empty((nq, nq, n_len))
	A_csum = numpy.empty((nq, nq, n_len))
	B = numpy.empty((t_max + 1, n_len))

	_p_value_backgrounds(f, A, B, A_csum, nq, n_bins, t_max, 
		numpy.uint64(offset))

	n = n_bins * nq + nq * offset
	for t in range(1, t_max + 1):
		assert numpy.all(B[t, :n] >= 0)
		assert numpy.all(B[t, :n] <= 1)
		assert numpy.all(numpy.diff(B[t, :n]) <= 1e-12)
		assert B[t, 0] > 0.9

	# Each single-column distribution sums to one.
	for i in range(nq):
		assert_array_almost_equal([A[i, i].sum()], [1.0], 12)


##


def test_tomtom_torch():
	# Tensors are converted with `.numpy()`, so results match numpy inputs
	# exactly, including float32 tensors against float32 arrays.
	torch = pytest.importorskip("torch")

	pwms = list(read_meme("tests/data/test.meme").values())
	pwms_t = [torch.from_numpy(pwm) for pwm in pwms]
	_assert_identical(tomtom(pwms, pwms), tomtom(pwms_t, pwms_t))
	_assert_identical(tomtom(pwms, pwms, n_nearest=3), 
		tomtom(pwms_t, pwms_t, n_nearest=3))

	pwms32 = [pwm.astype('float32') for pwm in pwms]
	_assert_identical(tomtom(pwms32, pwms), 
		tomtom([torch.from_numpy(pwm) for pwm in pwms32], pwms_t))


##


def _issue7_target(hi, lo):
	# The single target of issue #7: one near-one-hot column per position.
	t = numpy.full((4, 20), lo)
	for j, c in enumerate("ACCTACTAGGGGCTGAACCC"):
		t["ACGT".index(c), j] = hi

	return t


@pytest.mark.parametrize("hi, lo", [(0.997, 0.001), (0.9997, 0.0001)])
@pytest.mark.parametrize("reverse_complement", [True, False])
def test_tomtom_uniform_query_single_target(hi, lo, reverse_complement):
	# Every column of this target is the same distance from a uniform query
	# column, so every alignment scores the same and the p-value is 1. With
	# 0.997 the distances are exactly equal and the binned median divided by
	# a zero range. The second call checks that no error was left pending.
	q = numpy.full((4, 20), 0.25)
	T = _issue7_target(hi, lo)

	for _ in range(2):
		p = tomtom([q], [T], reverse_complement=reverse_complement)[0]
		assert_array_equal(p, [[1.0]])


@pytest.mark.parametrize("w", [6, 20])
def test_tomtom_uniform_query_round_off(w):
	# The same for 200 values of the dominant entry. When every column's
	# median is its minimum the scale was chosen from the largest shifted
	# score alone, which is round-off here, and p-values from 0.11 to 1 came
	# out, or the scale divided by zero.
	q = numpy.full((4, w), 0.25)

	for hi in numpy.linspace(0.5, 0.999, 200):
		T = _issue7_target(hi, (1 - hi) / 3)
		assert tomtom([q], [T])[0][0, 0] == 1


def test_tomtom_uniform_column_once():
	# A uniform column that occurs once in the query is scored inside the
	# parallel loop, where the division by zero returned uninitialized values
	# on the first call in a process and raised SystemError on later calls.
	q = numpy.random.RandomState(0).dirichlet(numpy.ones(4), size=12).T
	q[:, 5] = 0.25
	T = _issue7_target(0.997, 0.001)

	p = tomtom([q], [T])[0]
	assert_array_equal(p, tomtom([q], [T])[0])
	assert_array_almost_equal(p, [[0.8449]], 4)

	# A column 1e-9 from uniform never had a zero range, and scores the same.
	q[:, 5] = [0.25 + 1e-9, 0.25 - 1e-9, 0.25, 0.25]
	assert_array_almost_equal(p, tomtom([q], [T])[0], 4)


def test_tomtom_palindromic_single_column_target():
	# A one-column target that is its own reverse complement leaves a single
	# unique target column, so each query column has one distance and every
	# alignment scores the same.
	q = numpy.random.RandomState(1).dirichlet(numpy.ones(4), size=10).T
	T = numpy.array([[0.4], [0.1], [0.1], [0.4]])

	assert_array_equal(tomtom([q], [T])[0], [[1.0]])


@pytest.mark.parametrize("value", [-0.862561302169301, 0.0])
def test_binned_median_zero_range(value):
	# All values equal: the median is that value, with no division by the
	# zero range.
	X = numpy.full(7, value)
	counts = numpy.arange(1, 8, dtype='int64')
	bins = numpy.zeros((1000, 2), dtype='float64')
	zb = numpy.empty(7, dtype=numpy.int32)

	assert _binned_median(X, bins, value, value, counts) == value
	assert _binned_median_z(X, bins, value, value, counts, zb, 
		counts.sum() / 2) == value


def test_binned_median_block4_zero_range():
	# A constant row among three varied ones gets its value, and every row
	# matches `_binned_median_z` on its own.
	state = numpy.random.RandomState(2)
	n = 130
	rows = [numpy.full(n, -0.5)] + [state.uniform(-1.4, 0, n) for _ in range(3)]
	counts = state.randint(1, 5, size=n).astype('int64')
	halfway = counts.sum() / 2
	mn, mx = [r.min() for r in rows], [r.max() for r in rows]

	m = _binned_median_block4(rows[0], rows[1], rows[2], rows[3], mn[0], mx[0],
		mn[1], mx[1], mn[2], mx[2], mn[3], mx[3], numpy.zeros((1000, 2)), 
		counts, counts.astype(numpy.int32), numpy.empty((4, n), 
		dtype=numpy.int32), halfway)

	assert m[0] == -0.5
	for r, lo, hi, value in zip(rows, mn, mx, m):
		assert value == _binned_median_z(r, numpy.zeros((1000, 2)), lo, hi, 
			counts, numpy.empty(n, dtype=numpy.int32), halfway)


def test_p_value_backgrounds_lowest_bin():
	# With an offset of 0, a column's lowest score falls in bin 0, and that
	# bin's probability is part of the background: here S is 0 with
	# probability 0.25 and 3 with probability 0.75.
	f = numpy.zeros((1, 21))
	f[0, 0], f[0, 3] = 0.25, 0.75

	A = numpy.empty((1, 1, 40))
	A_csum = numpy.empty((1, 1, 40))
	B = numpy.empty((2, 40))
	_p_value_backgrounds(f, A, B, A_csum, 1, 20, 1, numpy.uint64(0))

	# B[1, j] is P(S >= j + 1).
	assert_array_almost_equal(B[1, :4], [0.75, 0.75, 0.75, 0.0])


def test_tomtom_one_hot_matches_meme():
	# One-hot queries against one-hot targets put most of each column's
	# probability in bin 0, which the background used to drop, and every
	# p-value came out as 1. The expected values are MEME 5.5.9's tomtom
	# (-dist ed -motif-pseudo 0) on the same motifs.
	Qs = one_hot_pwms(5, 5, 15, 30)[:3]
	Ts = one_hot_pwms(10, 5, 15, 31)

	p = tomtom(Qs, Ts)[0]
	assert_array_almost_equal(p, [
		[0.97555, 0.938792, 0.924108, 0.51685, 0.826821, 0.708033, 0.924108,
			0.708033, 0.897044, 0.215466],
		[0.075622, 0.041661, 0.84858, 0.787208, 0.649353, 0.899403, 0.84858,
			0.899403, 0.222739, 0.968174],
		[0.859372, 0.984423, 0.199181, 0.984423, 0.953196, 0.021359, 0.756237,
			0.997009, 0.548497, 0.398691]], 4)


def test_tomtom_self_matches_positive():
	# Each motif's match to itself lies in the far right tail of the
	# background, below the round-off of a cumsum that reaches 1, and
	# 1 - cumsum(pdf) gave 0 for every one of them.
	pwms = list(read_meme("tests/data/test.meme").values())
	p = tomtom(pwms, pwms)[0]

	assert (p > 0).all()
	assert (numpy.diag(p) < 1e-8).all()


def test_merge_rc_results_small_p():
	# 1 - (1 - p) ** 2 is 0 in float64 for p below about 1e-16.
	results = numpy.array([[1e-20, 5, 0, 5, 0], [1e-10, 3, 0, 5, 0],
		[1.0, 1, 0, 5, 0], [0.5, 2, 0, 5, 0]])
	_merge_rc_results(results)

	assert_allclose(results[:2, 0], [2e-20, 1e-10 * (2 - 1e-10)], rtol=1e-15)


def test_p_value_backgrounds_right_tail():
	# P(S >= s) is 1e-20 for s in [2, 19], which 1 - cumsum(pdf) returned as
	# 0 because the cumsum had already reached 1.
	f = numpy.zeros((1, 21))
	f[0, 1], f[0, 19] = 1.0, 1e-20

	A = numpy.empty((1, 1, 40))
	A_csum = numpy.empty((1, 1, 40))
	B = numpy.empty((2, 40))
	_p_value_backgrounds(f, A, B, A_csum, 1, 20, 1, numpy.uint64(0))

	# B[1, j] is P(S >= j + 1) for j < 20.
	assert B[1, 0] == 1
	assert_allclose(B[1, 1:19], 1e-20, rtol=1e-12)
	assert B[1, 19] == 0
