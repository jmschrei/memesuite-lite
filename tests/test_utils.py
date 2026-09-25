# test_utils.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import numpy
import pytest

from memelite.utils import characters
from memelite.utils import one_hot_encode
from memelite.utils import _fast_one_hot_encode

from numpy.testing import assert_raises
from numpy.testing import assert_array_almost_equal


##


def test_characters_ohe():
	seq = 'GCTAC'
	ohe = numpy.array([
		[0, 0, 0, 1, 0],
		[0, 1, 0, 0, 1],
		[1, 0, 0, 0, 0],
		[0, 0, 1, 0, 0]
	])

	seq_chars = characters(ohe)

	assert isinstance(seq_chars, str)
	assert len(seq_chars) == 5
	assert seq_chars == seq


def test_character_pwm():
	seq = 'GCTAC'
	ohe = numpy.array([
		[0.25, 0.00, 0.10, 0.95, 0.00],
		[0.20, 1.00, 1.00, 0.05, 1.00],
		[0.30, 0.00, 0.30, 0.00, 0.00],
		[0.25, 0.00, 3.00, 0.00, 0.00]
	])

	seq_chars = characters(ohe)

	assert isinstance(seq_chars, str)
	assert len(seq_chars) == 5
	assert seq_chars == seq


def test_characters_alphabet():
	seq = 'GCTAC'
	ohe = numpy.array([
		[0.25, 0.00, 0.10, 0.95, 0.00],
		[0.20, 1.00, 1.00, 0.05, 1.00],
		[0.30, 0.00, 0.30, 0.00, 0.00],
		[0.25, 0.00, 3.00, 0.00, 0.00],
		[0.00, 0.00, 0.00, 0.00, 0.00]
	])

	seq_chars = characters(ohe, ['A', 'C', 'G', 'T', 'N'])

	assert isinstance(seq_chars, str)
	assert len(seq_chars) == 5
	assert seq_chars == seq


def test_characters_raise_alphabet():
	seq = 'GCTAC'
	ohe = numpy.array([
		[0.25, 0.00, 0.10, 0.95, 0.00],
		[0.20, 1.00, 1.00, 0.05, 1.00],
		[0.30, 0.00, 0.30, 0.00, 0.00],
		[0.25, 0.00, 3.00, 0.00, 0.00]
	])

	assert_raises(ValueError, characters, ohe, ['A', 'C', 'G'])
	assert_raises(ValueError, characters, ohe, ['A', 'C', 'G', 'T', 'N'])


def test_characters_raise_dimensions():
	seq = 'GCTAC'
	#this will work for shape (1,4,5) but not for (N,4,5) where N > 1
	ohe = numpy.array([[
		[0.25, 0.00, 0.10, 0.95, 0.00],
		[0.20, 1.00, 1.00, 0.05, 1.00],
		[0.30, 0.00, 0.30, 0.00, 0.00],
		[0.25, 0.00, 3.00, 0.00, 0.00]
	]])
	
	assert characters(ohe) == seq
	
	ohe = numpy.concatenate([ohe, ohe], axis=0)
	assert_raises(ValueError, characters, ohe, ['A', 'C', 'G', 'T'])
	
	ohe = numpy.array([0.25, 0.00, 0.10, 0.95, 0.00])
	assert_raises(ValueError, characters, ohe, ['A', 'C', 'G', 'T'])


def test_characters_raise_ties():
	seq = 'GCTAC'
	ohe = numpy.array([
		[0.25, 0.00, 0.10, 0.95, 0.00],
		[0.20, 1.00, 1.00, 0.05, 1.00],
		[0.30, 1.00, 0.30, 0.00, 0.00],
		[0.25, 0.00, 3.00, 0.00, 0.00]
	])

	assert_raises(ValueError, characters, ohe, ['A', 'C', 'G', 'T'])
	assert characters(ohe, force=True) == seq


def test_characters_allow_N_all_zero():
	# An all-zero column should map to 'N' when allow_N is True.
	ohe = numpy.array([
		[1, 0, 0, 0],
		[0, 0, 0, 1],
		[0, 0, 1, 0],
		[0, 0, 0, 0]
	])

	assert characters(ohe, allow_N=True) == 'ANGC'


def test_characters_allow_N_all_columns_zero():
	# A PWM that is entirely zeros should be all 'N'.
	ohe = numpy.zeros((4, 5))

	assert characters(ohe, allow_N=True) == 'NNNNN'


def test_characters_allow_N_no_zero_columns():
	# When allow_N is True but no column is all-zero, the result matches
	# the standard argmax decoding.
	seq = 'GCTAC'
	ohe = numpy.array([
		[0, 0, 0, 1, 0],
		[0, 1, 0, 0, 1],
		[1, 0, 0, 0, 0],
		[0, 0, 1, 0, 0]
	])

	assert characters(ohe, allow_N=True) == seq


def test_characters_allow_N_resolves_ties():
	# allow_N bypasses the tie check entirely, so a tied column does not
	# raise and is resolved by argmax (earliest alphabet character).
	ohe = numpy.array([
		[1, 0, 0],
		[1, 0, 0],
		[0, 1, 0],
		[0, 0, 1]
	])

	assert characters(ohe, allow_N=True) == 'AGT'


def test_characters_allow_N_and_force():
	# With both force and allow_N set, all-zero columns become 'N' while
	# tied columns resolve to the earliest alphabet character.
	ohe = numpy.array([
		[1, 0, 0],
		[1, 0, 0],
		[0, 0, 1],
		[0, 0, 0]
	])

	assert characters(ohe, force=True, allow_N=True) == 'ANG'


def test_characters_force_tie_earliest_character():
	# When forcing, every tied position resolves to the character that
	# appears earliest in the alphabet.
	ohe = numpy.array([
		[1, 0],
		[1, 0],
		[0, 1],
		[0, 1]
	])

	assert characters(ohe, force=True) == 'AG'


def test_characters_3d_squeeze():
	# A (1, 4, L) input is squeezed to (4, L) and decoded normally.
	seq = 'GCTAC'
	ohe = numpy.array([[
		[0, 0, 0, 1, 0],
		[0, 1, 0, 0, 1],
		[1, 0, 0, 0, 0],
		[0, 0, 1, 0, 0]
	]])

	assert ohe.shape == (1, 4, 5)
	assert characters(ohe) == seq


def test_characters_single_column():
	# A single-column PWM returns a single character.
	ohe = numpy.array([[0], [1], [0], [0]])

	assert characters(ohe) == 'C'


def test_characters_custom_alphabet_ordering():
	# A reversed alphabet maps indices to the custom order.
	ohe = numpy.array([
		[1, 0],
		[0, 1],
		[0, 0],
		[0, 0]
	])

	assert characters(ohe, alphabet=['T', 'G', 'C', 'A']) == 'TG'


def test_characters_pwm_argmax():
	# Characterization test on a random PWM using float values; the decoded
	# sequence is the per-column argmax over the alphabet.
	random_state = numpy.random.RandomState(0)
	pwm = random_state.rand(4, 8)

	assert characters(pwm) == 'CGCTGCGA'


##


def test_one_hot_encode():
	seq = 'ACGTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 1],
		[0, 1, 0, 0, 0],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0]
	])
	seq_ohe = one_hot_encode(seq)

	assert seq_ohe.dtype == numpy.int8
	assert seq_ohe.shape == (4, 5)
	assert numpy.all(seq_ohe == ohe)

	seq = 'CCGTC'
	ohe = numpy.array([
		[0, 0, 0, 0, 0],
		[1, 1, 0, 0, 1],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0]
	])
	seq_ohe = one_hot_encode(seq)

	assert seq_ohe.dtype == numpy.int8
	assert seq_ohe.shape == (4, 5)
	assert numpy.all(seq_ohe == ohe)


	seq = 'AAAAA'
	ohe = numpy.array([
		[1, 1, 1, 1, 1],
		[0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0]
	])
	seq_ohe = one_hot_encode(seq)

	assert seq_ohe.dtype == numpy.int8
	assert seq_ohe.shape == (4, 5)
	assert numpy.all(seq_ohe == ohe)


def test_one_hot_encode_N():
	seq = 'ACGNNTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 0, 0, 1],
		[0, 1, 0, 0, 0, 0, 0],
		[0, 0, 1, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 1, 0]
	])
	seq_ohe = one_hot_encode(seq)

	assert seq_ohe.dtype == numpy.int8
	assert seq_ohe.shape == (4, 7)
	assert numpy.all(seq_ohe == ohe)

	seq = 'NNNNN'
	ohe = numpy.array([
		[0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0]
	])
	seq_ohe = one_hot_encode(seq)

	assert seq_ohe.dtype == numpy.int8
	assert seq_ohe.shape == (4, 5)
	assert numpy.all(seq_ohe == ohe)


def test_one_hot_encode_dtype():
	seq = 'ACGTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 1],
		[0, 1, 0, 0, 0],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0]
	])
	seq_ohe = one_hot_encode(seq, dtype=numpy.float32)

	assert seq_ohe.dtype == numpy.float32
	assert seq_ohe.shape == (4, 5)
	assert numpy.all(seq_ohe == ohe)


def test_one_hot_encode_alphabet():
	seq = 'ACGTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 1],
		[0, 1, 0, 0, 0],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0],
		[0, 0, 0, 0, 0]
	])
	seq_ohe = one_hot_encode(seq, alphabet=['A', 'C', 'G', 'T', 'Z'])
	
	assert seq_ohe.dtype == numpy.int8
	assert seq_ohe.shape == (5, 5)
	assert numpy.all(seq_ohe == ohe)


def test_one_hot_encode_ignore():
	seq = 'ACGTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 1],
		[0, 1, 0, 0, 0],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0],
		[0, 0, 0, 0, 0]
	])
	seq_ohe = one_hot_encode(seq, alphabet=['A', 'C', 'G', 'T', 'N'], ignore=[])
	
	assert seq_ohe.dtype == numpy.int8
	assert seq_ohe.shape == (5, 5)
	assert numpy.all(seq_ohe == ohe)


def test_one_hot_encode_alphabet():
	seq = 'ACGTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 1],
		[0, 1, 0, 0, 0],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0],
		[0, 0, 0, 0, 0]
	])
	seq_ohe = one_hot_encode(seq, alphabet=['A', 'C', 'G', 'T', 'Z'])
	
	assert seq_ohe.dtype == numpy.int8
	assert seq_ohe.shape == (5, 5)
	assert numpy.all(seq_ohe == ohe)


def test_one_hot_encode_raises_alphabet():
	seq = 'ACGTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 1],
		[0, 1, 0, 0, 0],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0]
	])

	assert_raises(ValueError, one_hot_encode, seq, ['A', 'C', 'G'])


def test_one_hot_encode_raises_ignore():
	seq = 'ACGTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 1],
		[0, 1, 0, 0, 0],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0],
		[0, 0, 0, 0, 0]
	])

	assert_raises(ValueError, one_hot_encode, seq, ['A', 'C', 'G', 'T', 'N'])


def test_one_hot_encode_lower_raises():
	seq = 'AcgTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 1],
		[0, 1, 0, 0, 0],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0]
	])

	assert_raises(ValueError, one_hot_encode, seq)


def test_one_hot_encode_unknown_char_raises():
	# A character that is neither in the alphabet nor in the ignore list
	# triggers a ValueError from the numba inner loop.
	assert_raises(ValueError, one_hot_encode, 'ACGTQ')
	assert_raises(ValueError, one_hot_encode, 'QACGT')
	assert_raises(ValueError, one_hot_encode, 'ACQGT')


def test_one_hot_encode_multi_char_ignore():
	# Multiple ignored characters all produce all-zero columns.
	seq = 'ACGTNX'
	ohe = numpy.array([
		[1, 0, 0, 0, 0, 0],
		[0, 1, 0, 0, 0, 0],
		[0, 0, 1, 0, 0, 0],
		[0, 0, 0, 1, 0, 0]
	])
	seq_ohe = one_hot_encode(seq, alphabet=['A', 'C', 'G', 'T'],
		ignore=['N', 'X'])

	assert seq_ohe.dtype == numpy.int8
	assert seq_ohe.shape == (4, 6)
	assert numpy.all(seq_ohe == ohe)
	assert numpy.all(seq_ohe.sum(axis=0) == numpy.array([1, 1, 1, 1, 0, 0]))


def test_one_hot_encode_ignore_all_zero_columns():
	# Ignored characters yield columns that sum to zero across the alphabet.
	seq = 'NANBN'
	seq_ohe = one_hot_encode(seq, alphabet=['A', 'B'], ignore=['N'])

	assert seq_ohe.shape == (2, 5)
	assert numpy.all(seq_ohe.sum(axis=0) == numpy.array([0, 1, 0, 1, 0]))


def test_one_hot_encode_dtype_float64():
	seq = 'ACGTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 1],
		[0, 1, 0, 0, 0],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0]
	])
	seq_ohe = one_hot_encode(seq, dtype=numpy.float64)

	assert seq_ohe.dtype == numpy.float64
	assert seq_ohe.shape == (4, 5)
	assert numpy.all(seq_ohe == ohe)


def test_one_hot_encode_dtype_bool():
	seq = 'ACGTA'
	ohe = numpy.array([
		[1, 0, 0, 0, 1],
		[0, 1, 0, 0, 0],
		[0, 0, 1, 0, 0],
		[0, 0, 0, 1, 0]
	], dtype=bool)
	seq_ohe = one_hot_encode(seq, dtype=bool)

	assert seq_ohe.dtype == numpy.bool_
	assert seq_ohe.shape == (4, 5)
	assert numpy.all(seq_ohe == ohe)


def test_one_hot_encode_empty_string():
	# An empty sequence produces an empty (alphabet_size, 0) encoding.
	seq_ohe = one_hot_encode('')

	assert seq_ohe.dtype == numpy.int8
	assert seq_ohe.shape == (4, 0)


def test_one_hot_encode_round_trip():
	# Decoding a one-hot encoding recovers the original sequence.
	seq = 'ACGTACGTGGCCAATTACGT'
	assert characters(one_hot_encode(seq)) == seq


def test_one_hot_encode_realistic_sequence():
	# Characterization test on a longer random sequence; the column sums are
	# all one (no ignored/unknown characters) and a round-trip recovers it.
	random_state = numpy.random.RandomState(0)
	seq = ''.join(numpy.array(['A', 'C', 'G', 'T'])[
		random_state.randint(0, 4, size=50)])

	seq_ohe = one_hot_encode(seq)

	assert seq_ohe.shape == (4, 50)
	assert numpy.all(seq_ohe.sum(axis=0) == 1)
	assert characters(seq_ohe) == seq
	assert seq[:10] == 'ATCATTTTCT'


def test_one_hot_encode_list_input():
	# one_hot_encode requires a str (it calls `bytearray(sequence, 'utf8')`),
	# as documented. Passing a list raises a TypeError.
	with pytest.raises(TypeError):
		one_hot_encode(['A', 'C', 'G', 'T'])


def test_fast_one_hot_encode_direct():
	# Direct test of the numba inner function: it fills the (n, m) matrix in
	# place using the byte->index mapping, skipping ignored indices (-1).
	seq = numpy.frombuffer(bytearray('ANC', 'utf8'), dtype=numpy.int8)
	mapping = numpy.zeros(256, dtype=numpy.int8) - 2
	mapping[ord('A')] = 0
	mapping[ord('C')] = 1
	mapping[ord('N')] = -1

	ohe = numpy.zeros((3, 2), dtype=numpy.int8)
	_fast_one_hot_encode(ohe, seq, mapping)

	expected = numpy.array([
		[1, 0],
		[0, 0],
		[0, 1]
	], dtype=numpy.int8)
	assert numpy.all(ohe == expected)


###

def _reference_one_hot(sequence, alphabet, ignore):
	"""A pure-Python reference encoding with shape (len(alphabet), len)."""

	ohe = numpy.zeros((len(alphabet), len(sequence)), dtype='int64')
	for i, char in enumerate(sequence):
		if char in ignore:
			continue

		ohe[list(alphabet).index(char), i] = 1

	return ohe


def _random_sequence(length, symbols, random_state):
	state = numpy.random.RandomState(random_state)
	return ''.join(numpy.array(list(symbols))[state.randint(0, len(symbols),
		size=length)])


_PROTEIN = list('ACDEFGHIKLMNPQRSTVWY')


@pytest.mark.parametrize("length", [0, 1, 2, 7, 100, 10000])
@pytest.mark.parametrize("alphabet,ignore", [
	(['A', 'C', 'G', 'T'], ['N']),
	(['A', 'B'], ['N']),
	(_PROTEIN, ['X', '-']),
	(['*', '-', '.', '0', '9'], ['?']),
	(['A', 'C', 'G', 'T'], [])
])
def test_one_hot_encode_reference(length, alphabet, ignore):
	symbols = list(alphabet) + list(ignore)
	seq = _random_sequence(length, symbols, random_state=length)

	ohe = one_hot_encode(seq, alphabet=alphabet, ignore=ignore)
	ohe_ref = _reference_one_hot(seq, alphabet, ignore)

	assert ohe.dtype == numpy.int8
	assert ohe.shape == (len(alphabet), length)
	assert numpy.array_equal(ohe, ohe_ref)

	# Every position sums to one except those holding an ignored character.
	is_ignored = numpy.array([c in ignore for c in seq], dtype=bool)
	assert numpy.all(ohe.sum(axis=0) == (~is_ignored).astype(int))


@pytest.mark.parametrize("length", [1, 2, 7, 100, 10000])
@pytest.mark.parametrize("alphabet", [['A', 'C', 'G', 'T'], ['A', 'B'],
	_PROTEIN, ['*', '-', '.', '0', '9']])
def test_one_hot_encode_characters_round_trip(length, alphabet):
	seq = _random_sequence(length, alphabet, random_state=length+1)
	ohe = one_hot_encode(seq, alphabet=alphabet, ignore=[])

	assert characters(ohe, alphabet=alphabet) == seq


@pytest.mark.parametrize("random_state", range(5))
def test_one_hot_encode_reference_with_N_round_trip(random_state):
	seq = _random_sequence(500, 'ACGTN', random_state=random_state)
	ohe = one_hot_encode(seq)

	assert numpy.array_equal(ohe, _reference_one_hot(seq, 'ACGT', 'N'))
	assert characters(ohe, allow_N=True) == seq


@pytest.mark.parametrize("dtype", [numpy.int8, numpy.uint8, numpy.int32,
	numpy.int64, numpy.float16, numpy.float32, numpy.float64, bool, 'float32',
	'int16'])
def test_one_hot_encode_dtypes(dtype):
	seq = _random_sequence(64, 'ACGTN', random_state=0)
	ohe = one_hot_encode(seq, dtype=dtype)

	assert ohe.dtype == numpy.dtype(dtype)
	assert ohe.shape == (4, 64)
	assert numpy.array_equal(ohe.astype('int64'),
		_reference_one_hot(seq, 'ACGT', 'N'))


def test_one_hot_encode_string_alphabet():
	# A string alphabet is treated like the equivalent list of characters.
	seq = _random_sequence(50, 'ACGT', random_state=0)
	assert numpy.array_equal(one_hot_encode(seq, alphabet='ACGT'),
		one_hot_encode(seq))


@pytest.mark.skip(reason="BUG: one_hot_encode documents tuple alphabets but "
	"raises TypeError because only lists are joined into a string")
def test_one_hot_encode_tuple_alphabet():
	seq = _random_sequence(50, 'ACGT', random_state=0)
	assert numpy.array_equal(one_hot_encode(seq, alphabet=('A', 'C', 'G', 'T')),
		one_hot_encode(seq))


def test_one_hot_encode_does_not_mutate_arguments():
	alphabet, ignore = ['A', 'C', 'G', 'T'], ['N']
	one_hot_encode('ACGTN', alphabet=alphabet, ignore=ignore)

	assert alphabet == ['A', 'C', 'G', 'T']
	assert ignore == ['N']


@pytest.mark.parametrize("position", [0, 1, 25, 49])
@pytest.mark.parametrize("char", ['Q', 'a', 'n', ' ', '\n', 'Z'])
def test_one_hot_encode_unknown_char_any_position(position, char):
	seq = list(_random_sequence(50, 'ACGT', random_state=position))
	seq[position] = char
	assert_raises(ValueError, one_hot_encode, ''.join(seq))


def test_one_hot_encode_non_ascii_raises():
	assert_raises(ValueError, one_hot_encode, 'ACGTé')


@pytest.mark.parametrize("alphabet,ignore", [
	(['A', 'C', 'G', 'T'], ['A']),
	(['A', 'C', 'G', 'T'], ['N', 'T']),
	(['A', 'C', 'G', 'T', 'N'], ['N']),
])
def test_one_hot_encode_raises_ignore_overlap(alphabet, ignore):
	assert_raises(ValueError, one_hot_encode, 'ACGT', alphabet=alphabet,
		ignore=ignore)


def test_one_hot_encode_empty_custom_alphabet():
	ohe = one_hot_encode('', alphabet=['A', 'B'])
	assert ohe.shape == (2, 0)


##


@pytest.mark.parametrize("random_state", range(10))
@pytest.mark.parametrize("n_chars", [2, 4, 20])
def test_characters_random_pwm_argmax(random_state, n_chars):
	# Continuous random values have no ties, so the result is the argmax.
	state = numpy.random.RandomState(random_state)
	alphabet = _PROTEIN[:n_chars]
	pwm = state.randn(n_chars, 37)

	expected = ''.join(alphabet[i] for i in pwm.argmax(axis=0))
	assert characters(pwm, alphabet=alphabet) == expected
	assert characters(pwm[None], alphabet=alphabet) == expected
	assert characters(pwm, alphabet=alphabet, force=True) == expected


def test_characters_negative_values():
	# Attribution-style matrices with only negative values still decode by
	# argmax.
	pwm = numpy.array([
		[-1.0, -2.0, -3.0],
		[-0.5, -3.0, -0.2],
		[-2.0, -0.1, -4.0],
		[-4.0, -5.0, -9.0]
	])

	assert characters(pwm) == 'CGC'
	assert characters(pwm, allow_N=True) == 'CGC'


@pytest.mark.parametrize("dtype", ['int8', 'int32', 'int64', 'uint8',
	'float16', 'float32', 'float64', 'bool'])
def test_characters_dtypes(dtype):
	seq = _random_sequence(40, 'ACGT', random_state=2)
	ohe = one_hot_encode(seq).astype(dtype)

	assert characters(ohe) == seq


def test_characters_tuple_alphabet():
	seq = _random_sequence(40, 'ACGT', random_state=3)
	assert characters(one_hot_encode(seq), alphabet=('A', 'C', 'G', 'T')) == seq


@pytest.mark.parametrize("random_state", range(5))
def test_characters_ties_force_allow_N(random_state):
	# Integer-valued PWMs with frequent ties and all-zero columns.
	state = numpy.random.RandomState(random_state)
	pwm = state.randint(0, 2, size=(4, 30))
	pwm[:, :3] = 0

	alphabet = numpy.array(['A', 'C', 'G', 'T'])
	is_max = pwm == pwm.max(axis=0, keepdims=True)
	assert is_max.sum(axis=0).max() > 1

	assert_raises(ValueError, characters, pwm)

	forced = ''.join(alphabet[pwm.argmax(axis=0)])
	assert characters(pwm, force=True) == forced

	with_N = ''.join('N' if pwm[:, i].sum() == 0 else c
		for i, c in enumerate(forced))
	assert characters(pwm, allow_N=True) == with_N
	assert characters(pwm, force=True, allow_N=True) == with_N
	assert with_N[:3] == 'NNN'


@pytest.mark.parametrize("shape", [(4,), (2, 4, 5), (1, 1, 4, 5), (3, 5),
	(5, 5), (1, 3, 5)])
def test_characters_raise_shapes(shape):
	pwm = numpy.random.RandomState(0).randn(*shape)
	assert_raises(ValueError, characters, pwm)
