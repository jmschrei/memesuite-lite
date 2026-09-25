# conftest.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import numba
import numpy
import pytest

from memelite.fimo import fimo
from memelite.tomtom import tomtom
from memelite.symmetric_tomtom import symmetric_tomtom
from memelite.utils import one_hot_encode


@pytest.fixture(scope="session", autouse=True)
def _jit_warmup():
	"""Compile the numba specializations the suite uses before any test runs.

	numba compiles one version of each kernel per combination of argument
	dtypes and memory layouts, which takes several seconds on a cold cache (a
	fresh clone, or any edit to a kernel). Without this, that cost lands on
	whichever test first uses a given combination. Warming up here keeps each
	test's own runtime a measure of the test and costs nothing extra overall.

	Each call below covers one specialization the suite reaches. Fortran-
	ordered PWMs (transposes, as `read_meme` returns) and C-ordered PWMs
	compile separately, as do hashed (`n_target_bins`) and unhashed targets,
	because hashing changes the dtype of the column index and the layout of
	the target matrix.
	"""

	# tomtom allocates a scratchpad per thread, so on a many-core machine the
	# default (every core) costs tens of GB and slows each call. Eight threads
	# still exercise the parallel paths.
	numba.set_num_threads(min(8, numba.config.NUMBA_NUM_THREADS))

	state = numpy.random.RandomState(0)
	pwms_f = [state.dirichlet(numpy.ones(4), size=n).T for n in (1, 4, 7)]
	pwms_c = [numpy.ascontiguousarray(pwm) for pwm in pwms_f]
	pwms_32 = [pwm.astype('float32') for pwm in pwms_f]
	onehot = [one_hot_encode("ACGTAC"), one_hot_encode("GGATC")]

	tomtom(pwms_f, pwms_f)
	tomtom(pwms_f, pwms_f, n_target_bins=None)
	tomtom(pwms_c, pwms_c)
	tomtom(pwms_c, pwms_c, n_target_bins=None)
	tomtom(pwms_32, pwms_32)
	tomtom(onehot, pwms_f)

	symmetric_tomtom(pwms_f)
	symmetric_tomtom(pwms_f, n_target_bins=None)
	symmetric_tomtom(pwms_c)
	symmetric_tomtom(pwms_32)

	X = numpy.eye(4, dtype='int8')[state.randint(4, size=(2, 30))].transpose(
		0, 2, 1)

	fimo({'a': pwms_f[2], 'b': pwms_f[1]}, "tests/data/test.fa")
	fimo({'a': pwms_f[2], 'b': pwms_f[1]}, X, bin_size=1)
	fimo({'a': pwms_c[2], 'b': pwms_c[1]}, X)
	fimo({'a': pwms_32[2], 'b': pwms_32[1]}, X)
