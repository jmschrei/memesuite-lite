# test_golden.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

"""Golden-output regression tests for tomtom, symmetric_tomtom and fimo.

The files in tests/data/golden/ pin the current outputs of the three public
functions over a grid of shapes, column contents and keyword arguments, which
is described by the `*_CASES` lists in tests/_golden_inputs.py. They exist so
that a rewrite of the numba kernels, e.g. during automated optimization, cannot
silently change any result.

Integer-valued outputs must match exactly: tomtom scores, offsets, overlaps,
strands and nearest-neighbor indexes, and the full fimo hit table (which motif,
sequence, start, end and strand). Floating point outputs are compared with a
relative tolerance of 1e-6.

A golden file pins behaviour, not correctness. Regenerate them with
`uv run --frozen python tests/generate_golden.py --force` only after an
intentional, reviewed change in behaviour, and say so in the commit message.
"""

import os
import numba
import numpy
import pytest

from memelite.tomtom import tomtom
from memelite.symmetric_tomtom import symmetric_tomtom
from memelite.fimo import fimo

from ._golden_inputs import GOLDEN_DIR
from ._golden_inputs import TOMTOM_CASES
from ._golden_inputs import SYMMETRIC_CASES
from ._golden_inputs import FIMO_CASES
from ._golden_inputs import build_tomtom
from ._golden_inputs import build_symmetric
from ._golden_inputs import build_fimo
from ._golden_inputs import fimo_table

from numpy.testing import assert_array_equal
from numpy.testing import assert_allclose


TOMTOM_KEYS = ['p', 'scores', 'offsets', 'overlaps', 'strands', 'idxs']
SYMMETRIC_KEYS = ['p', 'scores', 'offsets', 'overlaps', 'strands']

# Tomtom p-values of self-matches are round-off around zero (~1e-16), where a
# relative tolerance is meaningless, so they also get an absolute tolerance.
P_RTOL, P_ATOL = 1e-6, 1e-12

# fimo is compiled with fastmath, but each motif's scores are summed in a fixed
# order inside one prange iteration, so outputs are bitwise identical across
# runs and thread counts. The tolerance only absorbs formatting of the float
# columns and is far below the 0.1 score bin, so any change in the arithmetic
# that moves a score across a bin shows up as a changed p-value or hit set.
FIMO_RTOL = 1e-6


def _load(name):
	return numpy.load(os.path.join(GOLDEN_DIR, 'golden_{}.npz'.format(name)))


@pytest.fixture(scope='module')
def golden_tomtom():
	return _load('tomtom')


@pytest.fixture(scope='module')
def golden_symmetric():
	return _load('symmetric_tomtom')


@pytest.fixture(scope='module')
def golden_fimo():
	return _load('fimo')


def _case_ids(cases):
	return [case[0] for case in cases]


##


def test_golden_case_names_unique():
	for cases in (TOMTOM_CASES, SYMMETRIC_CASES, FIMO_CASES):
		names = _case_ids(cases)
		assert len(names) == len(set(names))


def test_golden_files_cover_cases(golden_tomtom, golden_symmetric,
	golden_fimo):
	# Every case has stored outputs and nothing stale is left in the files.
	for golden, cases in ((golden_tomtom, TOMTOM_CASES),
		(golden_symmetric, SYMMETRIC_CASES), (golden_fimo, FIMO_CASES)):
		stored = set(key.split('/')[0] for key in golden.files)
		assert stored == set(_case_ids(cases))


##


@pytest.mark.parametrize('n_jobs', [1, -1])
@pytest.mark.parametrize('case', TOMTOM_CASES, ids=_case_ids(TOMTOM_CASES))
def test_golden_tomtom(case, n_jobs, golden_tomtom):
	name = case[0]
	Qs, Ts, kwargs = build_tomtom(case)

	n_threads = numba.get_num_threads()
	out = tomtom(Qs, Ts, n_jobs=n_jobs, **kwargs)
	assert numba.get_num_threads() == n_threads

	n_outputs = 5 if kwargs.get('n_nearest') is None else 6
	assert len(out) == n_outputs

	for key, value in zip(TOMTOM_KEYS, out):
		expected = golden_tomtom['{}/{}'.format(name, key)]

		assert value.shape == expected.shape, key
		assert value.dtype == numpy.float64, key

		if key == 'p':
			assert_allclose(value, expected, rtol=P_RTOL, atol=P_ATOL,
				err_msg=key)
		else:
			assert_array_equal(value, expected, err_msg=key)


@pytest.mark.parametrize('n_jobs', [1, -1])
@pytest.mark.parametrize('case', SYMMETRIC_CASES,
	ids=_case_ids(SYMMETRIC_CASES))
def test_golden_symmetric_tomtom(case, n_jobs, golden_symmetric):
	name = case[0]
	Xs, kwargs = build_symmetric(case)

	out = symmetric_tomtom(Xs, n_jobs=n_jobs, **kwargs)
	assert len(out) == 5

	for key, value in zip(SYMMETRIC_KEYS, out):
		expected = golden_symmetric['{}/{}'.format(name, key)]
		value = value.copy()

		assert value.shape == expected.shape, key
		assert value.dtype == numpy.float64, key

		# The diagonal of offsets and overlaps is never written and holds
		# scratchpad garbage; the generator stores it as zero.
		if key in ('offsets', 'overlaps'):
			numpy.fill_diagonal(value, 0)

		if key == 'p':
			assert_allclose(value, expected, rtol=P_RTOL, atol=P_ATOL,
				err_msg=key)
		else:
			assert_array_equal(value, expected, err_msg=key)


@pytest.mark.parametrize('case', FIMO_CASES, ids=_case_ids(FIMO_CASES))
def test_golden_fimo(case, golden_fimo):
	name = case[0]
	motifs, sequences, kwargs = build_fimo(case)

	out = fimo(motifs, sequences, **kwargs)

	if kwargs.get('return_counts', False):
		expected = golden_fimo['{}/counts'.format(name)]
		assert isinstance(out, numpy.ndarray)
		assert out.dtype == expected.dtype
		assert_array_equal(out, expected)
		return

	table = fimo_table(out, dim=kwargs.get('dim', 0))
	for key, value in table.items():
		expected = golden_fimo['{}/{}'.format(name, key)]
		assert value.shape == expected.shape, key

		if key in ('score', 'pvalue'):
			assert_allclose(value, expected, rtol=FIMO_RTOL, err_msg=key)
		else:
			assert_array_equal(value, expected, err_msg=key)
