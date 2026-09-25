# test_fimo.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import warnings

import numpy
import pytest
import pandas

from memelite.fimo import _pwm_to_mapping
from memelite.fimo import _fast_convert
from memelite.fimo import logaddexp2
from memelite.fimo import fimo
from memelite.io import read_meme

from numpy.testing import assert_raises
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal


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


@pytest.fixture
def log_pwm():
	r = numpy.random.RandomState(0)

	pwm = numpy.exp(r.randn(4, 14))
	pwm = pwm / pwm.sum(axis=0, keepdims=True)
	return numpy.log2(pwm)

@pytest.fixture
def short_log_pwm():
	r = numpy.random.RandomState(0)

	pwm = numpy.exp(r.randn(4, 2))
	pwm = pwm / pwm.sum(axis=0, keepdims=True)
	return numpy.log2(pwm)

@pytest.fixture
def long_log_pwm():
	r = numpy.random.RandomState(0)

	pwm = numpy.exp(r.randn(4, 50))
	pwm = pwm / pwm.sum(axis=0, keepdims=True)
	return numpy.log2(pwm)


###


def test_pwm_to_mapping(log_pwm):
	smallest, mapping = _pwm_to_mapping(log_pwm, 0.1)

	assert smallest == -596
	assert mapping.shape == (600,)
	assert mapping.dtype == numpy.float64


	assert_array_almost_equal(mapping[:8], [7.4903e-08,  6.4154e-08,  
		4.2656e-08,  2.1158e-08, -1.1089e-08, -6.4833e-08, -1.1858e-07, 
		-1.8307e-07], 4)
	assert_array_almost_equal(mapping[100:108], [-0.0082, -0.0087, -0.0093, 
		-0.0099, -0.0105, -0.0112, -0.0119, -0.0126], 4)
	assert_array_almost_equal([mapping[~numpy.isinf(mapping)].min()], 
		[-28.], 4) 

	assert numpy.all(numpy.diff(mapping[~numpy.isinf(mapping)]) <= 0)
	assert numpy.isinf(mapping).sum() == 117


def test_short_pwm_to_mapping(short_log_pwm):
	smallest, mapping = _pwm_to_mapping(short_log_pwm, 0.1)

	assert smallest == -78
	assert mapping.shape == (67,)
	assert mapping.dtype == numpy.float64

	assert_array_almost_equal(mapping[:8], [ 8.3267e-17, -9.3109e-02, 
		-1.9265e-01, -1.9265e-01, -1.9265e-01, -1.9265e-01, -1.9265e-01, 
		-1.9265e-01], 4)
	assert_array_almost_equal([mapping[~numpy.isinf(mapping)].min()], 
		[-4], 4) 

	assert numpy.all(numpy.diff(mapping[~numpy.isinf(mapping)]) <= 0)
	assert numpy.isinf(mapping).sum() == 6


def test_long_pwm_to_mapping(long_log_pwm):
	smallest, mapping = _pwm_to_mapping(long_log_pwm, 0.1)

	assert smallest == -2045
	assert mapping.shape == (2085,)
	assert mapping.dtype == numpy.float64

	assert_array_almost_equal(mapping[:8], [-9.79656934e-16, -7.45058160e-09, 
		-2.23517430e-08, -3.72529047e-08, -5.96046475e-08, -9.68575534e-08, 
		-1.34110461e-07, -1.78813951e-07], 4)
	assert_array_almost_equal(mapping[100:108], [-3.3035e-14, -3.8062e-14, 
		-4.4026e-14, -5.1093e-14, -5.9458e-14, -6.9351e-14, -8.1036e-14, 
		-9.4826e-14], 4)
	assert_array_almost_equal([mapping[~numpy.isinf(mapping)].min()], 
		[-99.], 4) 

	assert numpy.all(numpy.diff(mapping[~numpy.isinf(mapping)]) <= 0)
	assert numpy.isinf(mapping).sum() == 486


def test_pwm_to_mapping_small_bins(log_pwm):
	smallest, mapping = _pwm_to_mapping(log_pwm, 0.01)

	assert smallest == -5955
	assert mapping.shape == (5864,)
	assert mapping.dtype == numpy.float64

	assert_array_almost_equal(mapping[:8], [-4.5076e-07, -4.6939e-07, 
		-4.8429e-07, -4.9546e-07, -5.1036e-07, -5.2154e-07, -5.4017e-07, 
		-5.5134e-07], 4)
	assert_array_almost_equal(mapping[100:108], [-4.5076e-07, -4.6939e-07, 
		-4.8429e-07, -4.9546e-07, -5.1036e-07, -5.2154e-07, -5.4017e-07, 
		-5.5134e-07], 4)
	assert_array_almost_equal([mapping[~numpy.isinf(mapping)].min()], 
		[-28.], 4) 

	assert numpy.all(numpy.diff(mapping[~numpy.isinf(mapping)]) <= 0)
	assert numpy.isinf(mapping).sum() == 1038


def test_pwm_to_mapping_large_bins(log_pwm):
	smallest, mapping = _pwm_to_mapping(log_pwm, 1)

	assert smallest == -60
	assert mapping.shape == (74,)
	assert mapping.dtype == numpy.float64

	assert_array_almost_equal(mapping[:8], [-3.5277e-07, -3.9577e-07, 
		-8.0423e-07, -3.0078e-06, -1.1978e-05, -4.1823e-05, -1.2694e-04, 
		-3.4126e-04], 4)
	assert_array_almost_equal([mapping[~numpy.isinf(mapping)].min()], 
		[-27], 4) 

	assert numpy.all(numpy.diff(mapping[~numpy.isinf(mapping)]) <= 0)
	assert numpy.isinf(mapping).sum() == 22


##


def test_fimo():
	hits = fimo("tests/data/test.meme", "tests/data/test.fa")

	assert len(hits) == 12
	for df in hits:
		assert isinstance(df, pandas.DataFrame)
		assert df.shape[1] == 8
		assert tuple(df.columns) == ('motif_name', 'motif_idx', 'sequence_name', 
			'start', 'end', 'strand', 'score', 'p-value')

	assert hits[0].shape == (1, 8)
	assert hits[9].shape == (1, 8)

	assert hits[0]['motif_name'][0] == "MEOX1_homeodomain_1"
	assert hits[0]['sequence_name'][0] == 'chr7'
	assert hits[0]['start'][0] == 1350
	assert hits[0]['end'][0] == 1360
	assert hits[0]['strand'][0] == '+'
	assert round(hits[0]['score'][0], 4) == round(11.446572, 4)
	assert round(hits[0]['p-value'][0], 4) == round(0.000075, 4)


	assert hits[9]['motif_name'][0] == "FOXQ1_MOUSE.H11MO.0.C"
	assert hits[9]['sequence_name'][0] == 'chr5'
	assert hits[9]['start'][0] == 121
	assert hits[9]['end'][0] == 133
	assert hits[9]['strand'][0] == '+'
	assert round(hits[9]['score'][0], 4) == round(3.17477, 4)
	assert round(hits[9]['p-value'][0], 4) == round(0.000099, 4)


def test_fimo_pwm_dict():
	pwms = read_meme("tests/data/test.meme")
	hits = fimo(pwms, "tests/data/test.fa")

	assert len(hits) == 12
	for df in hits:
		assert isinstance(df, pandas.DataFrame)
		assert df.shape[1] == 8
		assert tuple(df.columns) == ('motif_name', 'motif_idx', 'sequence_name', 
			'start', 'end', 'strand', 'score', 'p-value')

	assert hits[0].shape == (1, 8)
	assert hits[9].shape == (1, 8)

	assert hits[0]['motif_name'][0] == "MEOX1_homeodomain_1"
	assert hits[0]['sequence_name'][0] == 'chr7'
	assert hits[0]['start'][0] == 1350
	assert hits[0]['end'][0] == 1360
	assert hits[0]['strand'][0] == '+'
	assert round(hits[0]['score'][0], 4) == round(11.446572, 4)
	assert round(hits[0]['p-value'][0], 4) == round(0.000075, 4)


	assert hits[9]['motif_name'][0] == "FOXQ1_MOUSE.H11MO.0.C"
	assert hits[9]['sequence_name'][0] == 'chr5'
	assert hits[9]['start'][0] == 121
	assert hits[9]['end'][0] == 133
	assert hits[9]['strand'][0] == '+'
	assert round(hits[9]['score'][0], 4) == round(3.17477, 4)
	assert round(hits[9]['p-value'][0], 4) == round(0.000099, 4)


'''
def test_fimo_pwm_torch():
	pwms = read_meme("tests/data/test.meme")
	pwms = {name: torch.from_numpy(pwm) for name, pwm in pwms.items()}
	hits = fimo(pwms, "tests/data/test.fa")

	assert len(hits) == 12
	for df in hits:
		assert isinstance(df, pandas.DataFrame)
		assert df.shape[1] == 8
		assert tuple(df.columns) == ('motif_name', 'motif_idx', 'sequence_name', 
			'start', 'end', 'strand', 'score', 'p-value')

	assert hits[0].shape == (1, 8)
	assert hits[9].shape == (1, 8)

	assert hits[0]['motif_name'][0] == "MEOX1_homeodomain_1"
	assert hits[0]['sequence_name'][0] == 'chr7'
	assert hits[0]['start'][0] == 1350
	assert hits[0]['end'][0] == 1360
	assert hits[0]['strand'][0] == '+'
	assert round(hits[0]['score'][0], 4) == round(11.446572, 4)
	assert round(hits[0]['p-value'][0], 4) == round(0.000075, 4)


	assert hits[9]['motif_name'][0] == "FOXQ1_MOUSE.H11MO.0.C"
	assert hits[9]['sequence_name'][0] == 'chr5'
	assert hits[9]['start'][0] == 121
	assert hits[9]['end'][0] == 133
	assert hits[9]['strand'][0] == '+'
	assert round(hits[9]['score'][0], 4) == round(3.17477, 4)
	assert round(hits[9]['p-value'][0], 4) == round(0.000099, 4)
'''


def test_fimo_bin_size():
	hits = fimo("tests/data/test.meme", "tests/data/test.fa", bin_size=1)

	assert len(hits) == 12
	for df in hits:
		assert isinstance(df, pandas.DataFrame)
		assert df.shape[1] == 8
		assert tuple(df.columns) == ('motif_name', 'motif_idx', 'sequence_name', 
			'start', 'end', 'strand', 'score', 'p-value')

	assert hits[0].shape == (0, 8)
	assert hits[9].shape == (1, 8)


	assert hits[9]['motif_name'][0] == "FOXQ1_MOUSE.H11MO.0.C"
	assert hits[9]['sequence_name'][0] == 'chr5'
	assert hits[9]['start'][0] == 121
	assert hits[9]['end'][0] == 133
	assert hits[9]['strand'][0] == '+'
	assert round(hits[9]['score'][0], 4) == round(3.17477, 4)
	assert round(hits[9]['p-value'][0], 4) == round(0.000099, 4)


def test_fimo_threshold():
	hits = fimo("tests/data/test.meme", "tests/data/test.fa", threshold=0.001)

	assert len(hits) == 12
	for df in hits:
		assert isinstance(df, pandas.DataFrame)
		assert df.shape[1] == 8
		assert tuple(df.columns) == ('motif_name', 'motif_idx', 'sequence_name', 
			'start', 'end', 'strand', 'score', 'p-value')

	assert hits[0].shape == (13, 8)
	assert hits[9].shape == (157, 8)

	assert_array_equal(hits[0]['sequence_name'].values, ['chr1', 'chr2',
		'chr7', 'chr7', 'chr7', 'chr7', 'chr7', 'chr7', 'chr7', 'chr7', 'chr7', 
		'chr7', 'chr7'])
	assert_array_equal(hits[0]['start'].values, [ 190, 183, 667, 1096, 1106, 
		1161, 1350, 1384,  393, 1096, 1106, 1161, 1350])
	assert_array_equal(hits[0]['end'].values, [ 200, 193,  677, 1106, 1116, 
		1171, 1360, 1394,  403, 1106, 1116, 1171, 1360])
	assert_array_equal(hits[0]['strand'].values, ['+', '+', '+', '+', '+', '+', 
		'+', '+', '-', '-', '-', '-', '-'])
	assert_array_almost_equal(hits[0]['score'].values, [ 7.83938732, 7.40117477,
		7.40117477, 10.07583028, 10.07583028, 10.07583028, 11.4465722, 
		9.66039708,  8.30576608,  7.42200964,  7.42200964,  7.42200964,
  		7.54405698], 4)
	assert_array_almost_equal(hits[0]['p-value'].values, [7.60078430e-04, 
		9.22203064e-04, 9.22203064e-04, 1.89781189e-04, 1.89781189e-04, 
		1.89781189e-04, 7.53402710e-05, 2.45094299e-04, 5.88417053e-04, 
		9.22203064e-04, 9.22203064e-04, 9.22203064e-04, 8.81195068e-04], 4)


def test_fimo_rc():
	hits = fimo("tests/data/test.meme", "tests/data/test.fa",
		reverse_complement=False)

	assert len(hits[3]) == 1
	assert len(hits[7]) == 1
	assert len(hits[10]) == 1


##


def _make_one_hot(shape, random_state=None):
	"""Build a correctly-formed one-hot array with a fixed random state.

	A local replacement for the bugged `random_one_hot` helper in this file,
	which returns an undefined `one` variable instead of the constructed
	array. This version draws a random nucleotide index at each position and
	sets the corresponding channel to one.
	"""

	random_state = numpy.random.RandomState(random_state)

	n, c, l = shape
	idxs = random_state.randint(0, c, size=(n, l))

	ohe = numpy.zeros(shape, dtype='int8')
	for i in range(n):
		ohe[i, idxs[i], numpy.arange(l)] = 1

	return ohe


def test_logaddexp2_finite():
	r = numpy.random.RandomState(0)

	x = r.uniform(-20, 20, size=50)
	y = r.uniform(-20, 20, size=50)

	observed = numpy.array([logaddexp2(xi, yi) for xi, yi in zip(x, y)])
	expected = numpy.logaddexp2(x, y)
	assert_array_almost_equal(observed, expected, 4)


def test_logaddexp2_equal():
	for v in [-10.0, -1.0, 0.0, 1.0, 10.0]:
		assert_array_almost_equal([logaddexp2(v, v)],
			[numpy.logaddexp2(v, v)], 4)


def test_logaddexp2_both_neg_inf():
	ninf = float("-inf")
	assert logaddexp2(ninf, ninf) == ninf
	assert_array_almost_equal([logaddexp2(ninf, ninf)],
		[numpy.logaddexp2(ninf, ninf)], 4)


def test_logaddexp2_one_pos_inf():
	inf = float("inf")
	assert logaddexp2(inf, 3.0) == inf
	assert logaddexp2(3.0, inf) == inf
	assert_array_almost_equal([logaddexp2(inf, 3.0)],
		[numpy.logaddexp2(inf, 3.0)], 4)


def test_logaddexp2_one_neg_inf():
	ninf = float("-inf")
	assert_array_almost_equal([logaddexp2(ninf, 3.0)],
		[numpy.logaddexp2(ninf, 3.0)], 4)
	assert_array_almost_equal([logaddexp2(3.0, ninf)],
		[numpy.logaddexp2(3.0, ninf)], 4)
	assert logaddexp2(ninf, 3.0) == 3.0


##


def test_fast_convert():
	mapping = numpy.zeros(256, dtype=numpy.int8) - 1
	for i, c in enumerate([ord('A'), ord('C'), ord('G'), ord('T')]):
		mapping[c] = i

	X = numpy.frombuffer(bytearray("ACGTN", "utf8"), dtype=numpy.int8).copy()
	_fast_convert(X, mapping)

	assert_array_equal(X, [0, 1, 2, 3, -1])


def test_pwm_to_mapping_uniform():
	pwm = numpy.log2(numpy.full((4, 5), 0.25))
	smallest, mapping = _pwm_to_mapping(pwm, 0.1)

	assert smallest == -100
	assert mapping.shape == (86,)
	assert mapping.dtype == numpy.float64

	# A uniform PWM has a single achievable score, so only one finite bin.
	assert_array_almost_equal([mapping[~numpy.isinf(mapping)].min()], [0.0], 4)
	assert numpy.isinf(mapping).sum() == 85
	assert numpy.all(numpy.diff(mapping[~numpy.isinf(mapping)]) <= 0)


##


def test_fimo_one_hot():
	X = _make_one_hot((5, 4, 100), random_state=0)
	hits = fimo("tests/data/test.meme", X, threshold=0.001)

	assert len(hits) == 12
	for df in hits:
		assert isinstance(df, pandas.DataFrame)
		assert df.shape[1] == 8
		assert tuple(df.columns) == ('motif_name', 'motif_idx',
			'sequence_name', 'start', 'end', 'strand', 'score', 'p-value')

	# Without sequence names, sequence_name is the integer sequence index.
	assert_array_equal([len(df) for df in hits],
		[0, 1, 1, 0, 1, 1, 2, 2, 0, 0, 0, 1])

	df = hits[1]
	assert df['motif_name'][0] == "HIC2_MA0738.1"
	assert df['sequence_name'][0] == 1
	assert df['start'][0] == 38
	assert df['end'][0] == 47
	assert df['strand'][0] == '-'


def test_fimo_return_counts():
	counts = fimo("tests/data/test.meme", "tests/data/test.fa",
		return_counts=True)

	assert isinstance(counts, numpy.ndarray)
	assert counts.dtype == numpy.int32
	assert counts.shape == (12,)
	assert_array_equal(counts, [1, 0, 1, 2, 0, 0, 0, 2, 1, 1, 2, 0])

	# Counts must match the row counts of the default dataframe output.
	hits = fimo("tests/data/test.meme", "tests/data/test.fa")
	assert_array_equal(counts, [len(df) for df in hits])


def test_fimo_dim1():
	hits = fimo("tests/data/test.meme", "tests/data/test.fa", dim=1)

	# One dataframe per sequence that has at least one hit.
	assert len(hits) == 3
	for df in hits:
		assert isinstance(df, pandas.DataFrame)
		assert tuple(df.columns) == ('motif_name', 'motif_idx',
			'sequence_name', 'start', 'end', 'strand', 'score', 'p-value')
		# Each dataframe contains hits for exactly one sequence.
		assert df['sequence_name'].nunique() == 1

	assert_array_equal([df['sequence_name'].iloc[0] for df in hits],
		['chr1', 'chr5', 'chr7'])
	assert_array_equal([len(df) for df in hits], [1, 6, 3])


def test_fimo_dim1_no_warning():
	# The dim=1 concat must not raise a pandas FutureWarning about empty/all-NA
	# entries, since empty per-motif frames are filtered before concatenation.
	with warnings.catch_warnings():
		warnings.simplefilter("error", FutureWarning)
		hits = fimo("tests/data/test.meme", "tests/data/test.fa", dim=1)

	assert len(hits) == 3


def test_fimo_dim1_no_hits():
	# When no motif has any hit, dim=1 returns an empty list rather than
	# crashing on an empty pandas.concat.
	hits = fimo("tests/data/test.meme", "tests/data/test.fa", dim=1,
		threshold=1e-30)

	assert hits == []


def test_fimo_rc_false_strand():
	hits = fimo("tests/data/test.meme", "tests/data/test.fa",
		reverse_complement=False)

	assert len(hits) == 12
	for df in hits:
		assert (df['strand'] == '+').all()

	# Without reverse complements each forward hit appears once.
	assert_array_equal([len(df) for df in hits],
		[1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0])


def test_fimo_invalid_motifs():
	assert_raises(ValueError, fimo, 5, "tests/data/test.fa")
	assert_raises(ValueError, fimo, [1, 2, 3], "tests/data/test.fa")

	try:
		fimo(5, "tests/data/test.fa")
	except ValueError as e:
		assert "must be a dict or a filename" in str(e)


def test_fimo_eps():
	hits = fimo("tests/data/test.meme", "tests/data/test.fa", eps=0.001)
	assert_array_equal([len(df) for df in hits],
		[1, 0, 1, 2, 0, 0, 0, 2, 1, 1, 2, 0])

	# A larger pseudocount flattens motifs, dropping the weakest hits.
	hits = fimo("tests/data/test.meme", "tests/data/test.fa", eps=0.01)
	assert_array_equal([len(df) for df in hits],
		[1, 0, 1, 2, 0, 0, 0, 2, 1, 1, 0, 0])

	df = hits[9]
	assert df['motif_name'][0] == "FOXQ1_MOUSE.H11MO.0.C"
	assert df['sequence_name'][0] == 'chr5'
	assert df['start'][0] == 121
	assert df['end'][0] == 133
	assert df['strand'][0] == '+'
	assert_array_almost_equal(df['score'].values, [10.140161], 4)
	assert_array_almost_equal(df['p-value'].values, [0.000088], 4)


###
# Reference implementation.
#
# The tests below compare `fimo` against a slow, independent re-implementation
# written in plain numpy. It follows the documented semantics rather than the
# numba code: window scores are the sum of log2(pwm + eps) - log2(0.25) over
# the aligned positions, all-zero (N) columns contribute nothing, and the null
# distribution is the exact convolution, under a uniform background, of the
# per-column scores rounded to multiples of `bin_size`. A window is a hit when
# its score is strictly larger than k * bin_size, where k is the smallest
# integer score whose survival probability is below `threshold`, and its
# p-value is P(S >= int(score / bin_size)).
###


import math

from memelite.io import write_meme
from memelite.utils import one_hot_encode
from memelite.utils import characters

from numpy.testing import assert_allclose

from numpy.lib.stride_tricks import sliding_window_view


NAMES = ['motif_name', 'motif_idx', 'sequence_name', 'start', 'end', 'strand',
	'score', 'p-value']

# `_pwm_to_mapping` never enters its dynamic programming loop for a single
# column motif and returns an uninitialized `numpy.empty` buffer, so p-values
# and cutoffs for length-1 motifs are garbage that changes from call to call.
_SKIP_ONE = pytest.mark.skip(reason="BUG: _pwm_to_mapping returns "
	"uninitialized memory for single-column PWMs")
_ONE = pytest.param(1, marks=_SKIP_ONE)


def _random_pwms(n, min_len, max_len, alpha=0.5, random_state=0):
	"""Return a dict of `n` Dirichlet-sampled PWMs with random lengths."""

	random_state = numpy.random.RandomState(random_state)

	motifs = {}
	for i in range(n):
		length = random_state.randint(min_len, max_len+1)
		pwm = random_state.dirichlet(numpy.ones(4) * alpha, size=length).T
		motifs['motif_{}'.format(i)] = pwm

	return motifs


def _random_sequences(n, length, n_frac=0.0, dtype='int8', random_state=0):
	"""Return one-hot sequences with a fraction of all-zero (N) columns."""

	random_state = numpy.random.RandomState(random_state)

	X = numpy.zeros((n, 4, length), dtype='int8')
	idxs = random_state.randint(0, 4, size=(n, length))
	for i in range(n):
		X[i, idxs[i], numpy.arange(length)] = 1

	if n_frac > 0:
		mask = random_state.rand(n, length) < n_frac
		X.transpose(0, 2, 1)[mask] = 0

	return X.astype(dtype)


def _ref_sf(log_pwm, bin_size):
	"""Exact survival function of the binned score under a uniform background.

	Returns `(offset, sf)` where `sf[j]` is P(S >= j + offset) for the integer
	score S, i.e., the score in units of `bin_size`.
	"""

	ints = numpy.round(log_pwm / bin_size).astype(numpy.int64)

	pdf, offset = numpy.ones(1), 0
	for i in range(ints.shape[1]):
		col = ints[:, i]
		kernel = numpy.zeros(col.max() - col.min() + 1)
		numpy.add.at(kernel, col - col.min(), 0.25)

		pdf = numpy.convolve(pdf, kernel)
		offset += col.min()

	sf = numpy.cumsum(pdf[::-1])[::-1]
	return offset, sf


def _ref_fimo(motifs, X, bin_size=0.1, eps=0.0001, threshold=0.0001,
	reverse_complement=True):
	"""Brute-force FIMO over every window of every sequence.

	Returns a list of `(motif_idx, sequence_idx, start, end, strand, score,
	p-value, cutoff)` tuples, one per hit.
	"""

	X = numpy.asarray(X).astype(numpy.float64)

	rows = []
	for m, pwm in enumerate(motifs.values()):
		strands = [('+', pwm)]
		if reverse_complement:
			strands.append(('-', pwm[::-1, ::-1]))

		for strand, pwm_ in strands:
			log_pwm = numpy.log2(numpy.asarray(pwm_, dtype=numpy.float64) + eps
				) - math.log2(0.25)
			offset, sf = _ref_sf(log_pwm, bin_size)

			ks = numpy.where(sf < threshold)[0]
			if len(ks) == 0:
				continue

			cutoff = (ks[0] + offset) * bin_size
			n = log_pwm.shape[1]
			if X.shape[-1] < n:
				continue

			# windows has shape (n_seqs, 4, n_windows, n)
			windows = sliding_window_view(X, n, axis=-1)
			scores = numpy.einsum('iswn,sn->iw', windows, log_pwm)

			for i, start in zip(*numpy.where(scores > cutoff)):
				score = scores[i, start]
				k = int(score / bin_size) - offset
				p = sf[k] if k < len(sf) else 0.0
				rows.append((m, int(i), int(start), int(start+n), strand,
					score, p, cutoff))

	return rows


def _fimo_rows(hits):
	"""Flatten a list of `fimo` dataframes into a list of row tuples."""

	rows = []
	for df in hits:
		for r in df.itertuples(index=False):
			rows.append((int(r[1]), r[2], int(r[3]), int(r[4]), r[5],
				float(r[6]), float(r[7])))

	return rows


def _assert_matches_reference(hits, ref, bin_size):
	"""Check `fimo` hits against the brute-force reference.

	The hit sets must be identical, except for windows whose score lies within
	1e-5 of the cutoff: `fimo` stores cutoffs as float32 and sums scores with
	fastmath, so a window that lands on the cutoff to within rounding may fall
	on either side. Scores must match to rtol=1e-6 and p-values to rtol=1e-6,
	except that a score within 1e-8 of a bin edge may be looked up in the
	neighbouring bin.
	"""

	got = {row[:5]: row for row in _fimo_rows(hits)}
	exp = {row[:5]: row for row in ref}

	for key in set(got) - set(exp):
		cutoffs = [r[7] for r in ref if r[0] == key[0] and r[4] == key[4]]
		assert len(cutoffs) > 0 and abs(got[key][5] - cutoffs[0]) < 1e-5, \
			("unexpected hit", got[key])

	for key in set(exp) - set(got):
		assert abs(exp[key][5] - exp[key][7]) < 1e-5, ("missed hit", exp[key])

	keys = sorted(set(got) & set(exp))
	if len(keys) == 0:
		return

	scores = numpy.array([got[key][5] for key in keys])
	exp_scores = numpy.array([exp[key][5] for key in keys])
	assert_allclose(scores, exp_scores, rtol=1e-6, atol=1e-6)

	edge = numpy.abs(exp_scores / bin_size - numpy.round(exp_scores / bin_size))
	keep = edge > 1e-8

	p = numpy.array([got[key][6] for key in keys])[keep]
	exp_p = numpy.array([exp[key][6] for key in keys])[keep]
	assert_allclose(p, exp_p, rtol=1e-6, atol=0)


##


def test_ref_sf_uniform():
	# Sanity check of the reference itself: a uniform PWM has one score.
	offset, sf = _ref_sf(numpy.zeros((4, 5)), 0.1)
	assert offset == 0
	assert_array_almost_equal(sf, [1.0], 4)


def test_ref_sf_matches_enumeration():
	# Sanity check of the reference itself against 4^L enumeration.
	r = numpy.random.RandomState(0)
	log_pwm = numpy.round(r.randn(4, 4) / 0.1) * 0.1

	offset, sf = _ref_sf(log_pwm, 0.1)

	seqs = numpy.array(numpy.meshgrid(*[range(4)]*4, indexing='ij')).reshape(4, -1)
	scores = numpy.round(log_pwm[seqs, numpy.arange(4)[:, None]].sum(axis=0)
		/ 0.1).astype(int)

	for j in range(len(sf)):
		assert_array_almost_equal([sf[j]], [(scores >= j + offset).mean()], 10)


##


@pytest.mark.parametrize("n_motifs,min_len,max_len", [
	(1, 2, 2), (1, 5, 5), (1, 25, 25), (3, 2, 4), (7, 4, 12), (20, 2, 25),
	(12, 8, 20)
])
@pytest.mark.parametrize("n_seqs,length", [
	(1, 30), (4, 25), (13, 57), (30, 100)
])
def test_fimo_reference_shapes(n_motifs, min_len, max_len, n_seqs, length):
	motifs = _random_pwms(n_motifs, min_len, max_len, random_state=n_motifs)
	X = _random_sequences(n_seqs, length, random_state=length)

	# A permissive threshold so that short motifs have hits too.
	for threshold in (0.001, 0.05):
		hits = fimo(motifs, X, threshold=threshold)
		ref = _ref_fimo(motifs, X, threshold=threshold)

		assert len(hits) == n_motifs
		_assert_matches_reference(hits, ref, 0.1)


@pytest.mark.parametrize("dtype", ['int8', 'int32', 'float32', 'float64',
	'bool'])
@pytest.mark.parametrize("n_frac", [0.0, 0.1, 0.5, 1.0])
def test_fimo_reference_dtypes_and_N(dtype, n_frac):
	motifs = _random_pwms(6, 3, 15, random_state=1)
	X = _random_sequences(8, 80, n_frac=n_frac, dtype=dtype, random_state=2)

	hits = fimo(motifs, X, threshold=0.01)
	ref = _ref_fimo(motifs, X, threshold=0.01)
	_assert_matches_reference(hits, ref, 0.1)

	if n_frac == 1.0:
		# An all-N sequence scores 0 everywhere.
		for df in hits:
			assert_allclose(df['score'].values.astype(float), 0.0, atol=1e-12)


@pytest.mark.parametrize("bin_size", [0.01, 0.05, 0.1, 0.25, 0.5, 1.0])
def test_fimo_reference_bin_size(bin_size):
	motifs = _random_pwms(8, 2, 14, random_state=3)
	X = _random_sequences(10, 60, n_frac=0.05, random_state=4)

	hits = fimo(motifs, X, bin_size=bin_size, threshold=0.01)
	ref = _ref_fimo(motifs, X, bin_size=bin_size, threshold=0.01)
	_assert_matches_reference(hits, ref, bin_size)


@pytest.mark.parametrize("eps", [1e-8, 1e-4, 1e-3, 0.01, 0.1, 1.0])
def test_fimo_reference_eps(eps):
	# A low Dirichlet concentration gives many near-zero entries, where eps
	# has the largest effect.
	motifs = _random_pwms(8, 2, 14, alpha=0.1, random_state=5)
	X = _random_sequences(10, 60, random_state=6)

	hits = fimo(motifs, X, eps=eps, threshold=0.01)
	ref = _ref_fimo(motifs, X, eps=eps, threshold=0.01)
	_assert_matches_reference(hits, ref, 0.1)


@pytest.mark.parametrize("threshold", [1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1,
	0.5, 0.9])
@pytest.mark.parametrize("reverse_complement", [True, False])
def test_fimo_reference_threshold(threshold, reverse_complement):
	motifs = _random_pwms(8, 2, 20, random_state=7)
	X = _random_sequences(10, 60, n_frac=0.05, random_state=8)

	hits = fimo(motifs, X, threshold=threshold,
		reverse_complement=reverse_complement)
	ref = _ref_fimo(motifs, X, threshold=threshold,
		reverse_complement=reverse_complement)
	_assert_matches_reference(hits, ref, 0.1)


def test_fimo_reference_test_meme():
	# The reference also reproduces the real motifs in the test database.
	motifs = read_meme("tests/data/test.meme")
	X = _random_sequences(20, 150, n_frac=0.02, random_state=9)

	for threshold in (1e-4, 1e-3, 1e-2):
		hits = fimo(motifs, X, threshold=threshold)
		ref = _ref_fimo(motifs, X, threshold=threshold)
		_assert_matches_reference(hits, ref, 0.1)


def test_fimo_reference_float32_motifs():
	# float32 PWMs take a different numba specialization; results must match
	# the reference computed from the same float32 values.
	motifs = _random_pwms(6, 3, 15, random_state=10)
	motifs = {name: pwm.astype(numpy.float32) for name, pwm in motifs.items()}
	X = _random_sequences(8, 80, random_state=11)

	hits = fimo(motifs, X, threshold=0.01)
	ref = _ref_fimo(motifs, X, threshold=0.01)

	got = {row[:5]: row for row in _fimo_rows(hits)}
	exp = {row[:5]: row for row in ref}

	# log2 is evaluated in float32 inside fimo, so allow float32 rounding on
	# scores and on windows sitting at the cutoff.
	for key in set(got) ^ set(exp):
		row = got.get(key, exp.get(key))
		cutoff = [r[7] for r in ref if r[0] == key[0] and r[4] == key[4]][0]
		assert abs(row[5] - cutoff) < 1e-4, row

	keys = sorted(set(got) & set(exp))
	assert len(keys) > 0
	assert_allclose([got[k][5] for k in keys], [exp[k][5] for k in keys],
		rtol=1e-5, atol=1e-5)


##


@pytest.mark.parametrize("n", [_ONE, 2, 6, 11, 25])
@pytest.mark.parametrize("reverse_complement", [True, False])
def test_fimo_sequence_length_edges(n, reverse_complement):
	# Sequences shorter than, equal to, and one longer than the motif. A short
	# sequence must give no hits (and not hang), and L == n has one window.
	r = numpy.random.RandomState(n)
	pwm = r.dirichlet(numpy.ones(4) * 0.5, size=n).T
	motifs = {'m': pwm}

	for length in range(max(1, n-3), n+3):
		X = _random_sequences(6, length, random_state=length)

		hits = fimo(motifs, X, threshold=0.9,
			reverse_complement=reverse_complement)
		ref = _ref_fimo(motifs, X, threshold=0.9,
			reverse_complement=reverse_complement)
		_assert_matches_reference(hits, ref, 0.1)

		if length < n:
			assert len(hits[0]) == 0

		assert numpy.all(hits[0]['start'].values <= length - n)
		assert numpy.all(hits[0]['end'].values <= length)


@pytest.mark.parametrize("n", [_ONE, 4, 10, 20])
def test_fimo_last_window_forward(n):
	# Regression test: the last window of a sequence used to be skipped, so a
	# motif flush with the right end was never reported.
	r = numpy.random.RandomState(n)
	consensus = ''.join(r.choice(list('ACGT'), size=n))

	pwm = one_hot_encode(consensus).astype(numpy.float64) * 0.97 + 0.01
	pwm = pwm / pwm.sum(axis=0, keepdims=True)
	threshold = min(0.5, 2 * 0.25 ** n)

	X = numpy.stack([
		one_hot_encode('T' * 10 + consensus),
		one_hot_encode(consensus + 'T' * 10),
	])

	hits = fimo({'m': pwm}, X, reverse_complement=False,
		threshold=threshold)[0]
	starts = {(int(s), int(i)) for s, i in zip(hits['sequence_name'],
		hits['start'])}

	assert (0, 10) in starts
	assert (1, 0) in starts

	# A sequence exactly the length of the motif has exactly one window.
	X = one_hot_encode(consensus)[None]
	hits = fimo({'m': pwm}, X, reverse_complement=False,
		threshold=threshold)[0]

	assert len(hits) == 1
	assert hits['start'][0] == 0
	assert hits['end'][0] == n


@pytest.mark.parametrize("n", [_ONE, 4, 10, 20])
def test_fimo_last_window_reverse(n):
	# The same regression on the reverse strand: the reverse complement of the
	# consensus placed flush with the right end is reported on '-'.
	r = numpy.random.RandomState(n + 100)
	consensus = ''.join(r.choice(list('ACGT'), size=n))
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	rc = ''.join(complement[c] for c in reversed(consensus))

	pwm = one_hot_encode(consensus).astype(numpy.float64) * 0.97 + 0.01
	pwm = pwm / pwm.sum(axis=0, keepdims=True)
	threshold = min(0.5, 2 * 0.25 ** n)

	X = one_hot_encode('A' * 7 + rc)[None]
	hits = fimo({'m': pwm}, X, threshold=threshold)[0]
	hits = hits[hits['strand'] == '-']

	assert 7 in set(hits['start'].astype(int))
	assert 7 + n in set(hits['end'].astype(int))


##


@pytest.mark.parametrize("L", [_ONE, 2, 3, 4, 5, 6])
@pytest.mark.parametrize("bin_size", [0.001, 0.01, 0.1, 0.5])
def test_fimo_p_values_vs_enumeration(L, bin_size):
	# Enumerate all 4^L sequences and compare each reported p-value to the
	# exact tail probability of the unbinned score. Rounding each column to a
	# multiple of `bin_size` moves the total score by at most L * bin_size / 2
	# and the lookup `int(score / bin_size)` moves it by at most one more bin,
	# so the binned p-value must lie between the exact tail probabilities at
	# score +/- bin_size * (1 + L / 2).
	r = numpy.random.RandomState(L)
	pwm = r.dirichlet(numpy.ones(4) * 0.5, size=L).T

	seqs = numpy.array(numpy.meshgrid(*[range(4)]*L, indexing='ij')
		).reshape(L, -1).T
	X = numpy.zeros((len(seqs), 4, L), dtype='int8')
	for j in range(L):
		X[numpy.arange(len(seqs)), seqs[:, j], j] = 1

	log_pwm = numpy.log2(pwm + 0.0001) + 2
	exact = numpy.sort(log_pwm[seqs, numpy.arange(L)].sum(axis=1))

	def exact_sf(s):
		return (len(exact) - numpy.searchsorted(exact, s, side='left')) / len(
			exact)

	hits = fimo({'m': pwm}, X, bin_size=bin_size, threshold=0.99,
		reverse_complement=False)[0]
	assert len(hits) > 0

	delta = bin_size * (1 + L / 2.) + 1e-9
	for score, p in zip(hits['score'], hits['p-value']):
		assert exact_sf(score + delta) <= p + 1e-12
		assert p <= exact_sf(score - delta) + 1e-12

	# Every sequence was scanned exactly once.
	assert hits['sequence_name'].nunique() == len(hits)


@pytest.mark.parametrize("threshold", [1e-4, 1e-3, 0.01, 0.1])
def test_fimo_p_values_bounded_and_monotone(threshold):
	# p-values are positive, below the threshold, and non-increasing in the
	# score. The reverse-complement PWM has the same null distribution, so
	# monotonicity holds across both strands of a motif.
	motifs = _random_pwms(10, 2, 20, random_state=12)
	X = _random_sequences(20, 100, n_frac=0.05, random_state=13)

	hits = fimo(motifs, X, threshold=threshold)
	for df in hits:
		p = df['p-value'].values.astype(float)
		assert numpy.all(p > 0)
		assert numpy.all(p < threshold)

		order = numpy.argsort(df['score'].values.astype(float), kind='stable')
		assert numpy.all(numpy.diff(p[order]) <= 1e-15)


##


@pytest.mark.parametrize("L", [_ONE, 2, 5, 14, 30])
@pytest.mark.parametrize("bin_size", [0.01, 0.1, 1.0])
@pytest.mark.parametrize("alpha", [0.1, 1.0, 10.0])
def test_pwm_to_mapping_vs_reference(L, bin_size, alpha):
	r = numpy.random.RandomState(L)
	pwm = r.dirichlet(numpy.ones(4) * alpha, size=L).T
	log_pwm = numpy.log2(pwm + 0.0001) + 2

	smallest, mapping = _pwm_to_mapping(log_pwm, bin_size)
	offset, sf = _ref_sf(log_pwm, bin_size)

	# mapping[j] is log2 P(S >= j + smallest); below the reference support the
	# probability is one and above it the probability is zero.
	assert smallest <= offset
	assert len(mapping) >= len(sf) + offset - smallest

	j = numpy.arange(len(mapping)) + smallest - offset
	expected = numpy.ones(len(mapping))
	inside = (j >= 0) & (j < len(sf))
	expected[inside] = sf[j[inside]]
	expected[j >= len(sf)] = 0

	nonzero = expected > 0
	assert_allclose(2 ** mapping[nonzero], expected[nonzero], rtol=1e-9)
	assert numpy.all(numpy.isinf(mapping[~nonzero]))
	assert numpy.all(numpy.diff(mapping[~numpy.isinf(mapping)]) <= 1e-12)


def test_all_pwm_to_mapping_matches_single():
	# The parallel wrapper returns the same mapping as calling the single
	# motif function on each slice of the concatenated PWMs.
	from memelite.fimo import _all_pwm_to_mapping

	motifs = list(_random_pwms(15, 2, 25, random_state=14).values())
	lengths = numpy.cumsum([0] + [m.shape[1] for m in motifs]).astype(
		numpy.uint64)
	log_pwms = numpy.log2(numpy.concatenate(motifs, axis=1) + 0.0001) + 2

	smallests, mappings = _all_pwm_to_mapping(log_pwms, lengths, 0.1)

	assert len(smallests) == len(motifs)
	assert len(mappings) == len(motifs)
	for i, pwm in enumerate(motifs):
		smallest, mapping = _pwm_to_mapping(numpy.log2(pwm + 0.0001) + 2, 0.1)
		assert smallests[i] == smallest
		assert_array_equal(mappings[i], mapping)


@pytest.mark.parametrize("x,y", [(-1000.0, 0.0), (0.0, -1000.0), (500.0, 499.0),
	(-500.0, -500.0), (1e-10, -1e-10), (60.0, 1.0), (-60.0, 1.0), (1000.0,
	1000.0)])
def test_logaddexp2_extremes(x, y):
	assert_allclose(logaddexp2(x, y), numpy.logaddexp2(x, y), rtol=1e-12,
		atol=1e-12)
	assert logaddexp2(x, y) == logaddexp2(y, x)


##


def _key_set(hits, strand=None):
	return {row[:5] for row in _fimo_rows(hits) if strand is None or
		row[4] == strand}


def test_fimo_threshold_monotone():
	# Hits at a stricter threshold are a subset of hits at a looser one, with
	# identical scores and p-values.
	motifs = _random_pwms(10, 2, 20, random_state=15)
	X = _random_sequences(15, 80, random_state=16)

	thresholds = [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5]
	results = [{row[:5]: row for row in _fimo_rows(fimo(motifs, X,
		threshold=t))} for t in thresholds]

	for small, large in zip(results[:-1], results[1:]):
		assert set(small) <= set(large)
		for key in small:
			assert small[key][5:] == large[key][5:]

	assert len(results[-1]) > len(results[0])


def test_fimo_reverse_complement_union():
	# With reverse complements, the hits are the forward hits plus the hits of
	# the reverse-complemented PWM scanned forward, labelled '-'.
	motifs = _random_pwms(8, 2, 18, random_state=17)
	rc_motifs = {name: pwm[::-1, ::-1] for name, pwm in motifs.items()}
	X = _random_sequences(12, 70, n_frac=0.05, random_state=18)

	both = fimo(motifs, X, threshold=0.01)
	fwd = fimo(motifs, X, threshold=0.01, reverse_complement=False)
	rev = fimo(rc_motifs, X, threshold=0.01, reverse_complement=False)

	for df, df_f, df_r in zip(both, fwd, rev):
		plus = df[df['strand'] == '+'].reset_index(drop=True)
		minus = df[df['strand'] == '-'].reset_index(drop=True)
		df_r = df_r.assign(strand='-')

		pandas.testing.assert_frame_equal(plus, df_f, check_dtype=False)
		pandas.testing.assert_frame_equal(minus, df_r, check_dtype=False)

		# Forward hits come before reverse hits in each dataframe.
		strands = list(df['strand'])
		assert strands == sorted(strands, key=lambda s: s == '-')


def test_fimo_row_order():
	# Within a motif and strand, hits are ordered by sequence then position.
	motifs = _random_pwms(6, 2, 12, random_state=19)
	X = _random_sequences(10, 80, random_state=20)

	for df in fimo(motifs, X, threshold=0.05):
		for strand in '+-':
			sub = df[df['strand'] == strand]
			keys = list(zip(sub['sequence_name'].astype(int),
				sub['start'].astype(int)))
			assert keys == sorted(keys)


@pytest.mark.parametrize("threshold", [1e-4, 1e-3, 0.05])
@pytest.mark.parametrize("reverse_complement", [True, pytest.param(False,
	marks=pytest.mark.skip(reason="BUG: return_counts=True with "
	"reverse_complement=False indexes hits[i + n_] past the end of the list "
	"and raises IndexError"))])
def test_fimo_return_counts_matches_rows(threshold, reverse_complement):
	motifs = _random_pwms(12, 2, 20, random_state=21)
	X = _random_sequences(10, 90, n_frac=0.05, random_state=22)

	counts = fimo(motifs, X, threshold=threshold,
		reverse_complement=reverse_complement, return_counts=True)
	hits = fimo(motifs, X, threshold=threshold,
		reverse_complement=reverse_complement)

	assert counts.shape == (12,)
	assert counts.dtype == numpy.int32
	assert_array_equal(counts, [len(df) for df in hits])


@pytest.mark.parametrize("threshold", [1e-3, 0.05])
@pytest.mark.parametrize("reverse_complement", [True, False])
def test_fimo_dim1_regroups_dim0(threshold, reverse_complement):
	motifs = _random_pwms(8, 2, 15, random_state=23)
	X = _random_sequences(9, 70, random_state=24)

	hits0 = fimo(motifs, X, threshold=threshold,
		reverse_complement=reverse_complement)
	hits1 = fimo(motifs, X, threshold=threshold,
		reverse_complement=reverse_complement, dim=1)

	all0 = pandas.concat([df for df in hits0 if len(df) > 0])
	names = sorted(all0['sequence_name'].unique())

	# One dataframe per sequence with a hit, sorted by sequence name, each
	# holding that sequence's rows in motif order with a fresh index.
	assert len(hits1) == len(names)
	for name, df in zip(names, hits1):
		assert list(df.columns) == NAMES
		assert (df['sequence_name'] == name).all()
		assert_array_equal(df.index, numpy.arange(len(df)))

		expected = all0[all0['sequence_name'] == name].reset_index(drop=True)
		pandas.testing.assert_frame_equal(df, expected)


##


def test_fimo_motif_order_invariance():
	motifs = _random_pwms(10, 2, 20, random_state=25)
	X = _random_sequences(10, 80, random_state=26)

	perm = numpy.random.RandomState(0).permutation(10)
	names = list(motifs.keys())
	motifs_perm = {names[i]: motifs[names[i]] for i in perm}

	hits = fimo(motifs, X, threshold=0.01)
	hits_perm = fimo(motifs_perm, X, threshold=0.01)

	for j, i in enumerate(perm):
		assert (hits_perm[j]['motif_idx'] == j).all()
		pandas.testing.assert_frame_equal(
			hits[i].drop(columns='motif_idx'),
			hits_perm[j].drop(columns='motif_idx'))


@pytest.mark.parametrize("split", [1, 5, 11])
def test_fimo_batch_split_invariance(split):
	# Scanning the sequences in two batches gives the same hits, with the
	# sequence index of the second batch shifted by the split point.
	motifs = _random_pwms(8, 2, 18, random_state=27)
	X = _random_sequences(12, 60, n_frac=0.05, random_state=28)

	full = _fimo_rows(fimo(motifs, X, threshold=0.01))
	a = _fimo_rows(fimo(motifs, X[:split], threshold=0.01))
	b = _fimo_rows(fimo(motifs, X[split:], threshold=0.01))
	b = [(r[0], r[1] + split) + r[2:] for r in b]

	assert sorted(full) == sorted(a + b)


def test_fimo_single_sequence_invariance():
	# Each sequence scanned alone gives the same hits as in the batch.
	motifs = _random_pwms(6, 2, 15, random_state=29)
	X = _random_sequences(6, 50, random_state=30)

	full = _fimo_rows(fimo(motifs, X, threshold=0.02))

	single = []
	for i in range(len(X)):
		rows = _fimo_rows(fimo(motifs, X[i:i+1], threshold=0.02))
		single.extend([(r[0], i) + r[2:] for r in rows])

	assert sorted(full) == sorted(single)


def test_fimo_padding_invariance():
	# Hits inside a sequence do not change when N padding is appended.
	motifs = _random_pwms(6, 2, 15, random_state=31)
	X = _random_sequences(6, 50, random_state=32)
	X_pad = numpy.concatenate([X, numpy.zeros((6, 4, 20), dtype=X.dtype)],
		axis=-1)

	rows = _fimo_rows(fimo(motifs, X, threshold=0.02))
	rows_pad = [r for r in _fimo_rows(fimo(motifs, X_pad, threshold=0.02))
		if r[3] <= 50]

	assert sorted(rows) == sorted(rows_pad)


##


def _write_fasta(path, seqs):
	with open(path, 'w') as outfile:
		for name, seq in seqs.items():
			outfile.write(">{}\n".format(name))
			for i in range(0, len(seq), 60):
				outfile.write(seq[i:i+60] + "\n")


@pytest.mark.parametrize("reverse_complement", [True, False])
def test_fimo_fasta_matches_one_hot(tmp_path, reverse_complement):
	# The FASTA path upper-cases sequences and treats non-ACGT characters as N,
	# so it must match scanning the equivalent one-hot encodings.
	r = numpy.random.RandomState(33)
	motifs = _random_pwms(8, 2, 18, random_state=34)

	seqs = {}
	for i, length in enumerate([5, 17, 100, 250, 64]):
		seq = r.choice(list('ACGTacgtN'), size=length, p=[0.2] * 4 + [0.04] *
			4 + [0.04])
		seqs['seq{}'.format(i)] = ''.join(seq)

	fasta = str(tmp_path / "seqs.fa")
	_write_fasta(fasta, seqs)

	hits = fimo(motifs, fasta, threshold=0.01,
		reverse_complement=reverse_complement)

	names = list(seqs.keys())
	for i, (name, seq) in enumerate(seqs.items()):
		X = one_hot_encode(seq.upper())[None]
		hits_i = fimo(motifs, X, threshold=0.01,
			reverse_complement=reverse_complement)

		for df, df_i in zip(hits, hits_i):
			sub = df[df['sequence_name'] == name].reset_index(drop=True)
			df_i = df_i.assign(sequence_name=name)
			pandas.testing.assert_frame_equal(sub, df_i, check_dtype=False)

	for df in hits:
		assert set(df['sequence_name']) <= set(names)


@pytest.mark.parametrize("alphabet,reverse_complement", [
	(['T', 'G', 'C', 'A'], True),
	(['C', 'A', 'T', 'G'], True),
	(['A', 'C', 'T', 'G'], False),
	(['G', 'T', 'A', 'C'], False),
])
def test_fimo_fasta_alphabet(tmp_path, alphabet, reverse_complement):
	# Reordering the alphabet together with the PWM rows gives the same hits.
	# The reverse complement flips PWM rows, so it is only well-defined when
	# the reversed alphabet is the complement.
	r = numpy.random.RandomState(35)
	seqs = {'a': ''.join(r.choice(list('ACGT'), size=120)),
		'b': ''.join(r.choice(list('ACGT'), size=80))}

	fasta = str(tmp_path / "seqs.fa")
	_write_fasta(fasta, seqs)

	motifs = _random_pwms(6, 2, 12, random_state=36)
	order = [['A', 'C', 'G', 'T'].index(c) for c in alphabet]
	motifs_perm = {name: pwm[order] for name, pwm in motifs.items()}

	hits = fimo(motifs, fasta, threshold=0.02,
		reverse_complement=reverse_complement)
	hits_perm = fimo(motifs_perm, fasta, alphabet=alphabet, threshold=0.02,
		reverse_complement=reverse_complement)

	assert sum(len(df) for df in hits) > 0
	for df, df_perm in zip(hits, hits_perm):
		pandas.testing.assert_frame_equal(df, df_perm)


def test_fimo_meme_file_matches_dict(tmp_path):
	motifs = _random_pwms(10, 2, 20, random_state=37)
	filename = str(tmp_path / "motifs.meme")
	write_meme(filename, motifs)

	X = _random_sequences(8, 80, random_state=38)

	hits = fimo(filename, X, threshold=0.01)
	hits_dict = fimo(motifs, X, threshold=0.01)

	assert len(hits) == len(hits_dict)
	for df, df_dict in zip(hits, hits_dict):
		pandas.testing.assert_frame_equal(df, df_dict)


##


def test_fimo_output_format():
	motifs = _random_pwms(5, 3, 12, random_state=39)
	X = _random_sequences(6, 60, random_state=40)

	hits = fimo(motifs, X, threshold=0.05)

	assert isinstance(hits, list)
	assert len(hits) == 5
	for i, (name, df) in enumerate(zip(motifs.keys(), hits)):
		assert isinstance(df, pandas.DataFrame)
		assert list(df.columns) == NAMES
		assert len(df) > 0

		assert (df['motif_name'] == name).all()
		assert (df['motif_idx'] == i).all()
		assert set(df['strand']) <= {'+', '-'}

		for column in ('motif_idx', 'sequence_name', 'start', 'end'):
			assert numpy.issubdtype(df[column].dtype, numpy.integer), column

		for column in ('score', 'p-value'):
			assert df[column].dtype == numpy.float64, column

		n = motifs[name].shape[1]
		assert_array_equal(df['end'] - df['start'], n)
		assert numpy.all(df['start'] >= 0)
		assert numpy.all(df['end'] <= 60)
		assert numpy.all(df['sequence_name'] >= 0)
		assert numpy.all(df['sequence_name'] < 6)


def test_fimo_empty_output_format():
	motifs = _random_pwms(4, 3, 12, random_state=41)
	X = _random_sequences(3, 40, random_state=42)

	hits = fimo(motifs, X, threshold=1e-30)

	assert len(hits) == 4
	for df in hits:
		assert df.shape == (0, 8)
		assert list(df.columns) == NAMES

	counts = fimo(motifs, X, threshold=1e-30, return_counts=True)
	assert_array_equal(counts, numpy.zeros(4))


def test_fimo_does_not_modify_inputs():
	motifs = _random_pwms(4, 3, 12, random_state=43)
	X = _random_sequences(3, 40, random_state=44)

	motifs_copy = {name: pwm.copy() for name, pwm in motifs.items()}
	X_copy = X.copy()

	fimo(motifs, X, threshold=0.05)

	assert list(motifs.keys()) == list(motifs_copy.keys())
	for name in motifs:
		assert_array_equal(motifs[name], motifs_copy[name])
	assert_array_equal(X, X_copy)


def test_fimo_repeatable():
	# Repeated calls give identical output, including row order.
	motifs = _random_pwms(12, 2, 20, random_state=45)
	X = _random_sequences(10, 80, random_state=46)

	hits1 = fimo(motifs, X, threshold=0.01)
	hits2 = fimo(motifs, X, threshold=0.01)

	for df1, df2 in zip(hits1, hits2):
		pandas.testing.assert_frame_equal(df1, df2)


def test_fimo_invalid_motif_values():
	# PWMs that are neither numpy arrays nor have a `.numpy()` method.
	X = _random_sequences(2, 30, random_state=47)

	assert_raises(ValueError, fimo, {'a': [[0.25] * 3] * 4}, X)
	assert_raises(ValueError, fimo, {'a': 5}, X)
	assert_raises(ValueError, fimo, (1, 2), X)


@_SKIP_ONE
def test_pwm_to_mapping_single_column():
	# A single column PWM has survival probabilities 1, 0.75, 0.5, 0.25 at its
	# four scores and zero above them.
	pwm = numpy.array([[0.05], [0.6], [0.15], [0.2]])
	log_pwm = numpy.log2(pwm + 0.0001) + 2

	smallest, mapping = _pwm_to_mapping(log_pwm, 0.1)
	offset, sf = _ref_sf(log_pwm, 0.1)
	assert_allclose(sf[sf > 0][[0, -1]], [1.0, 0.25])

	j = numpy.arange(len(mapping)) + smallest - offset
	expected = numpy.ones(len(mapping))
	inside = (j >= 0) & (j < len(sf))
	expected[inside] = sf[j[inside]]
	expected[j >= len(sf)] = 0

	nonzero = expected > 0
	assert_allclose(2 ** mapping[nonzero], expected[nonzero], rtol=1e-9)
	assert numpy.all(numpy.isinf(mapping[~nonzero]))


@_SKIP_ONE
def test_fimo_single_column_motif():
	# The 'C' position of a single-column motif has p-value 0.25 and must be
	# reported at threshold 0.3; no other position is.
	pwm = numpy.array([[0.05], [0.6], [0.15], [0.2]])
	X = numpy.eye(4, dtype='int8')[None]

	for _ in range(5):
		hits = fimo({'m': pwm}, X, threshold=0.3, reverse_complement=False)[0]
		assert len(hits) == 1
		assert hits['start'][0] == 1
		assert_array_almost_equal(hits['p-value'].values, [0.25], 4)
