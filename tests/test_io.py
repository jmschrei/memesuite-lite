# test_io.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import os

import numpy
import pytest

from memelite.io import read_meme
from memelite.io import write_meme

from numpy.testing import assert_raises
from numpy.testing import assert_array_almost_equal


TEST_MEME = os.path.join(os.path.dirname(__file__), "data", "test.meme")
TEST2_MEME = os.path.join(os.path.dirname(__file__), "data", "test2.meme")


TEST_NAMES = [
	'MEOX1_homeodomain_1',
	'HIC2_MA0738.1',
	'GCR_HUMAN.H11MO.0.A',
	'FOSL2+JUND_MA1145.1',
	'TEAD3_TEA_2',
	'ZN263_HUMAN.H11MO.0.A',
	'PAX7_PAX_2',
	'SMAD3_MA0795.1',
	'MEF2D_HUMAN.H11MO.0.A',
	'FOXQ1_MOUSE.H11MO.0.C',
	'TBX19_MA0804.1',
	'Hes1_MA1099.1'
]

TEST2_NAMES = [
	'MA0636.1 BHLHE41',
	'MA0641.1 ELF4',
	'MA0042.2 FOXI1',
	'MA0033.2 FOXL1'
]


### read_meme


def test_read_meme_returns_dict():
	motifs = read_meme(TEST_MEME)

	assert isinstance(motifs, dict)
	assert len(motifs) == 12


def test_read_meme_names():
	motifs = read_meme(TEST_MEME)
	assert list(motifs.keys()) == TEST_NAMES


def test_read_meme_shapes():
	motifs = read_meme(TEST_MEME)

	for name, pwm in motifs.items():
		assert pwm.ndim == 2
		assert pwm.shape[0] == 4

	assert motifs['MEOX1_homeodomain_1'].shape == (4, 10)
	assert motifs['HIC2_MA0738.1'].shape == (4, 9)


def test_read_meme_values():
	motifs = read_meme(TEST_MEME)
	pwm = motifs['MEOX1_homeodomain_1']

	# Columns are positions; rows are A, C, G, T.
	assert_array_almost_equal(pwm[:, 0],
		[0.30190, 0.23459, 0.32856, 0.13495], 4)
	assert_array_almost_equal(pwm[:, 4],
		[0.99834, 0.00000, 0.00067, 0.00100], 4)
	assert_array_almost_equal(pwm[:, 9],
		[0.19060, 0.43619, 0.20460, 0.16861], 4)

	# Columns of a PWM should sum to ~1.
	assert_array_almost_equal(pwm.sum(axis=0), numpy.ones(10), 4)

	pwm2 = motifs['HIC2_MA0738.1']
	assert_array_almost_equal(pwm2[:, 0],
		[0.46108, 0.06971, 0.40428, 0.06493], 4)


def test_read_meme_n_motifs_one():
	motifs = read_meme(TEST_MEME, n_motifs=1)

	assert len(motifs) == 1
	assert list(motifs.keys()) == ['MEOX1_homeodomain_1']
	assert motifs['MEOX1_homeodomain_1'].shape == (4, 10)


def test_read_meme_n_motifs_middle():
	motifs = read_meme(TEST_MEME, n_motifs=5)

	assert len(motifs) == 5
	assert list(motifs.keys()) == TEST_NAMES[:5]


def test_read_meme_n_motifs_none():
	motifs = read_meme(TEST_MEME, n_motifs=None)

	assert len(motifs) == 12
	assert list(motifs.keys()) == TEST_NAMES


def test_read_meme_n_motifs_too_large():
	motifs = read_meme(TEST_MEME, n_motifs=20)

	assert len(motifs) == 12
	assert list(motifs.keys()) == TEST_NAMES


def test_read_meme_test2():
	motifs = read_meme(TEST2_MEME)

	assert len(motifs) == 4
	assert list(motifs.keys()) == TEST2_NAMES

	assert motifs['MA0636.1 BHLHE41'].shape == (4, 10)
	assert motifs['MA0641.1 ELF4'].shape == (4, 12)

	assert_array_almost_equal(motifs['MA0636.1 BHLHE41'][:, 0],
		[0.302309, 0.109266, 0.582263, 0.006163], 4)
	assert_array_almost_equal(motifs['MA0636.1 BHLHE41'][:, 9],
		[0.008609, 0.648441, 0.114588, 0.228362], 4)


def test_read_meme_test2_n_motifs():
	assert len(read_meme(TEST2_MEME, n_motifs=1)) == 1
	assert len(read_meme(TEST2_MEME, n_motifs=2)) == 2
	assert len(read_meme(TEST2_MEME, n_motifs=4)) == 4
	assert len(read_meme(TEST2_MEME, n_motifs=10)) == 4
	assert len(read_meme(TEST2_MEME, n_motifs=None)) == 4


### write_meme


def test_write_meme_header(tmp_path):
	motifs = read_meme(TEST_MEME, n_motifs=1)
	filename = str(tmp_path / "out.meme")

	write_meme(filename, motifs)

	with open(filename, "r") as infile:
		contents = infile.read()

	assert "MEME version 4\n" in contents
	assert "ALPHABET= ACGT\n" in contents
	assert "strands: + -\n" in contents
	assert "Background letter frequencies\n" in contents
	assert "A 0.25 C 0.25 G 0.25 T 0.25\n" in contents


def test_write_meme_dict_names(tmp_path):
	motifs = read_meme(TEST_MEME)
	filename = str(tmp_path / "out.meme")

	write_meme(filename, motifs)

	with open(filename, "r") as infile:
		names = [line[6:].strip() for line in infile if line[:5] == 'MOTIF']

	assert names == TEST_NAMES


def test_write_meme_list_names(tmp_path):
	motifs = list(read_meme(TEST_MEME).values())
	filename = str(tmp_path / "out.meme")

	write_meme(filename, motifs)

	with open(filename, "r") as infile:
		names = [line[6:].strip() for line in infile if line[:5] == 'MOTIF']

	assert names == [str(i) for i in range(12)]


def test_write_meme_lpm_line(tmp_path):
	# The letter-probability matrix line should be written from pwm.shape,
	# so alength comes from shape[0] and w comes from shape[1].
	motifs = read_meme(TEST_MEME, n_motifs=2)
	filename = str(tmp_path / "out.meme")

	write_meme(filename, motifs)

	with open(filename, "r") as infile:
		lpm = [line.strip() for line in infile if line[:6] == 'letter']

	assert lpm[0] == \
		"letter-probability matrix: alength= 4 w= 10 nsites= 1 E= 0"
	assert lpm[1] == \
		"letter-probability matrix: alength= 4 w= 9 nsites= 1 E= 0"


### round-trip


def test_round_trip_dict(tmp_path):
	motifs = read_meme(TEST_MEME)
	filename = str(tmp_path / "out.meme")

	write_meme(filename, motifs)
	motifs2 = read_meme(filename)

	assert list(motifs.keys()) == list(motifs2.keys())

	# Python's default float formatting round-trips exactly, so the recovered
	# PWMs match the originals to full precision.
	for name in motifs:
		assert motifs[name].shape == motifs2[name].shape
		assert_array_almost_equal(motifs[name], motifs2[name], 4)


def test_round_trip_list(tmp_path):
	motifs = read_meme(TEST2_MEME)
	pwms = list(motifs.values())
	filename = str(tmp_path / "out.meme")

	write_meme(filename, pwms)
	motifs2 = read_meme(filename)

	assert list(motifs2.keys()) == [str(i) for i in range(len(pwms))]

	for pwm, pwm2 in zip(pwms, motifs2.values()):
		assert pwm.shape == pwm2.shape
		assert_array_almost_equal(pwm, pwm2, 4)


### edge cases


def test_single_motif_file(tmp_path):
	motifs = read_meme(TEST_MEME, n_motifs=1)
	filename = str(tmp_path / "single.meme")

	write_meme(filename, motifs)
	motifs2 = read_meme(filename)

	assert len(motifs2) == 1
	assert list(motifs2.keys()) == ['MEOX1_homeodomain_1']
	assert_array_almost_equal(
		motifs['MEOX1_homeodomain_1'],
		motifs2['MEOX1_homeodomain_1'], 4)


def test_differing_widths(tmp_path):
	# Construct motifs with different widths and confirm each round-trips
	# with its own width preserved.
	pwm_a = numpy.full((4, 6), 0.25)
	pwm_b = numpy.full((4, 11), 0.25)
	motifs = {"short": pwm_a, "long": pwm_b}

	filename = str(tmp_path / "widths.meme")
	write_meme(filename, motifs)
	motifs2 = read_meme(filename)

	assert motifs2["short"].shape == (4, 6)
	assert motifs2["long"].shape == (4, 11)
	assert_array_almost_equal(motifs2["short"], pwm_a, 4)
	assert_array_almost_equal(motifs2["long"], pwm_b, 4)


### randomized round-trips


def _random_pwms(n, min_len=1, max_len=50, random_state=0):
	state = numpy.random.RandomState(random_state)

	pwms = []
	for i in range(n):
		length = state.randint(min_len, max_len+1)
		pwm = state.dirichlet(numpy.ones(4) * 0.5, size=length).T
		pwms.append(pwm)

	return pwms


@pytest.mark.parametrize("n", [1, 2, 7, 100])
@pytest.mark.parametrize("random_state", [0, 1, 2])
def test_round_trip_random_dict_exact(tmp_path, n, random_state):
	pwms = _random_pwms(n, random_state=random_state)
	motifs = {"motif_{}".format(i): pwm for i, pwm in enumerate(pwms)}
	filename = str(tmp_path / "out.meme")

	write_meme(filename, motifs)
	motifs2 = read_meme(filename)

	assert list(motifs2.keys()) == list(motifs.keys())

	# `write_meme` formats floats with repr, which round-trips exactly.
	for name, pwm in motifs.items():
		assert motifs2[name].shape == pwm.shape
		assert motifs2[name].dtype == numpy.float64
		assert numpy.array_equal(motifs2[name], pwm)


@pytest.mark.parametrize("n", [1, 5, 50])
def test_round_trip_random_list_exact(tmp_path, n):
	pwms = _random_pwms(n, random_state=3)
	filename = str(tmp_path / "out.meme")

	write_meme(filename, pwms)
	motifs2 = read_meme(filename)

	assert list(motifs2.keys()) == [str(i) for i in range(n)]
	for pwm, pwm2 in zip(pwms, motifs2.values()):
		assert numpy.array_equal(pwm, pwm2)


@pytest.mark.parametrize("width", [1, 2, 3, 8, 21, 50])
def test_round_trip_widths(tmp_path, width):
	pwm = _random_pwms(1, min_len=width, max_len=width, random_state=width)[0]
	filename = str(tmp_path / "out.meme")

	write_meme(filename, {"m": pwm})
	motifs2 = read_meme(filename)

	assert motifs2["m"].shape == (4, width)
	assert numpy.array_equal(motifs2["m"], pwm)


def test_round_trip_extreme_values(tmp_path):
	# Exact zeros, exact ones, subnormal-scale and many-digit values.
	pwm = numpy.array([
		[0.0, 1.0, 0.1234567890123456, 1e-300, 0.25],
		[1.0, 0.0, 0.8765432109876544, 1 - 1e-16, 0.25],
		[0.0, 0.0, 0.0, 5e-324, 0.25],
		[0.0, 0.0, 0.0, 0.0, 0.25]
	])
	filename = str(tmp_path / "out.meme")

	write_meme(filename, {"m": pwm})
	motifs2 = read_meme(filename)

	assert numpy.array_equal(motifs2["m"], pwm)


@pytest.mark.parametrize("name", ["simple", "with space", "MA0001.1 AGL3",
	"a/b:c|d", "FOSL2+JUND", "x-y_z(1)", "weird[]{};'\""])
def test_round_trip_names(tmp_path, name):
	pwm = numpy.full((4, 3), 0.25)
	filename = str(tmp_path / "out.meme")

	write_meme(filename, {name: pwm})
	motifs2 = read_meme(filename)

	assert list(motifs2.keys()) == [name]


def test_round_trip_n_motifs_random(tmp_path):
	pwms = _random_pwms(30, random_state=4)
	filename = str(tmp_path / "out.meme")
	write_meme(filename, pwms)

	for k in [1, 2, 15, 29, 30, 31]:
		motifs2 = read_meme(filename, n_motifs=k)
		assert list(motifs2.keys()) == [str(i) for i in range(min(k, 30))]

		for pwm, pwm2 in zip(pwms, motifs2.values()):
			assert numpy.array_equal(pwm, pwm2)


### write_meme exact text


def test_write_meme_golden_text(tmp_path):
	pwm = numpy.array([[0.5, 0.0], [0.25, 1.0], [0.125, 0.0], [0.125, 0.0]])
	filename = str(tmp_path / "out.meme")

	write_meme(filename, {"t": pwm})

	with open(filename, "r") as infile:
		contents = infile.read()

	assert contents == ("MEME version 4\n\n"
		"ALPHABET= ACGT\n\n"
		"strands: + -\n\n"
		"Background letter frequencies\n"
		"A 0.25 C 0.25 G 0.25 T 0.25\n\n"
		"MOTIF t\n"
		"letter-probability matrix: alength= 4 w= 2 nsites= 1 E= 0\n"
		"0.5 0.25 0.125 0.125\n"
		"0.0 1.0 0.0 0.0\n"
		"URL BLANK\n\n")


def test_write_meme_rows_per_motif(tmp_path):
	pwms = _random_pwms(10, random_state=5)
	filename = str(tmp_path / "out.meme")
	write_meme(filename, pwms)

	with open(filename, "r") as infile:
		lines = infile.read().split("\n")

	starts = [i for i, line in enumerate(lines) if line.startswith("MOTIF")]
	assert len(starts) == 10

	for start, pwm in zip(starts, pwms):
		assert lines[start+1] == ("letter-probability matrix: alength= 4 "
			"w= {} nsites= 1 E= 0".format(pwm.shape[1]))

		rows = lines[start+2:start+2+pwm.shape[1]]
		assert all(len(row.split()) == 4 for row in rows)
		assert lines[start+2+pwm.shape[1]] == "URL BLANK"


def test_write_meme_empty(tmp_path):
	filename = str(tmp_path / "out.meme")
	write_meme(filename, {})

	assert read_meme(filename) == {}
	write_meme(filename, [])
	assert read_meme(filename) == {}


### read_meme format robustness


_HEADER = "MEME version 4\n\nALPHABET= ACGT\n\n"
_MOTIF_A = ("MOTIF a\n"
	"letter-probability matrix: alength= 4 w= 2 nsites= 1 E= 0\n"
	"0.1 0.2 0.3 0.4\n"
	"0.25 0.25 0.25 0.25\n")
_MOTIF_B = ("MOTIF b\n"
	"letter-probability matrix: alength= 4 w= 1 nsites= 1 E= 0\n"
	"1 0 0 0\n")

_PWM_A = numpy.array([[0.1, 0.25], [0.2, 0.25], [0.3, 0.25], [0.4, 0.25]])
_PWM_B = numpy.array([[1.0], [0.0], [0.0], [0.0]])


def _read_text(tmp_path, text, **kwargs):
	filename = str(tmp_path / "in.meme")
	with open(filename, "w", newline="") as outfile:
		outfile.write(text)

	return read_meme(filename, **kwargs)


def _assert_ab(motifs):
	assert list(motifs.keys()) == ['a', 'b']
	assert numpy.array_equal(motifs['a'], _PWM_A)
	assert numpy.array_equal(motifs['b'], _PWM_B)


def test_read_meme_standard_text(tmp_path):
	motifs = _read_text(tmp_path, _HEADER + _MOTIF_A + "URL x\n\n" + _MOTIF_B
		+ "URL y\n")
	_assert_ab(motifs)


def test_read_meme_crlf(tmp_path):
	text = _HEADER + _MOTIF_A + "URL x\n\n" + _MOTIF_B + "URL y\n"
	motifs = _read_text(tmp_path, text.replace("\n", "\r\n"))
	_assert_ab(motifs)


def test_read_meme_extra_blank_lines(tmp_path):
	motifs = _read_text(tmp_path, _HEADER + "\n\n" + _MOTIF_A + "\n\n\n"
		+ _MOTIF_B + "\n\n\n")
	_assert_ab(motifs)


def test_read_meme_tab_separated(tmp_path):
	text = _MOTIF_A.replace("0.1 0.2 0.3 0.4", "0.1\t0.2\t0.3\t0.4")
	motifs = _read_text(tmp_path, _HEADER + text + "\n")

	assert numpy.array_equal(motifs['a'], _PWM_A)


def test_read_meme_row_whitespace(tmp_path):
	text = _MOTIF_A.replace("0.1 0.2 0.3 0.4", "  0.1  0.2 0.3 0.4  ")
	motifs = _read_text(tmp_path, _HEADER + text + "\n")

	assert numpy.array_equal(motifs['a'], _PWM_A)


def test_read_meme_scientific_notation(tmp_path):
	text = _MOTIF_A.replace("0.1 0.2", "1e-1 2.0E-01")
	motifs = _read_text(tmp_path, _HEADER + text + "\n")

	assert numpy.array_equal(motifs['a'], _PWM_A)


def test_read_meme_lpm_extra_fields(tmp_path):
	text = _MOTIF_A.replace("E= 0", "E= 1.2e-5 extra 7")
	motifs = _read_text(tmp_path, _HEADER + text + "\n")

	assert numpy.array_equal(motifs['a'], _PWM_A)


def test_read_meme_alternate_name(tmp_path):
	# The key is the full remainder of the MOTIF line, including an alternate
	# name when one is given.
	text = _MOTIF_A.replace("MOTIF a", "MOTIF MA0001.1 AGL3")
	motifs = _read_text(tmp_path, _HEADER + text + "\n")

	assert list(motifs.keys()) == ['MA0001.1 AGL3']


@pytest.mark.parametrize("line", ["MOTIF\tMA0001.1\tAGL3",
	"MOTIF  MA0001.1  AGL3", "MOTIF MA0001.1 AGL3 \t", "MOTIF \tMA0001.1 \t AGL3",
	"MOTIFMA0001.1 AGL3"])
def test_read_meme_header_whitespace(tmp_path, line):
	# MEME allows any whitespace, or none, after MOTIF and any whitespace
	# between the fields. The key joins the fields with one space.
	text = _MOTIF_A.replace("MOTIF a", line)
	motifs = _read_text(tmp_path, _HEADER + text)

	assert list(motifs.keys()) == ['MA0001.1 AGL3']
	assert numpy.array_equal(motifs['MA0001.1 AGL3'], _PWM_A)


def test_read_meme_empty_file(tmp_path):
	assert _read_text(tmp_path, "") == {}


def test_read_meme_header_only(tmp_path):
	assert _read_text(tmp_path, _HEADER) == {}


@pytest.mark.parametrize("trailing", ["", "\n"])
def test_read_meme_ends_after_matrix(tmp_path, trailing):
	text = _HEADER + _MOTIF_A + "\n" + _MOTIF_B.rstrip("\n") + trailing
	_assert_ab(_read_text(tmp_path, text))


def test_read_meme_no_separator(tmp_path):
	text = _HEADER + _MOTIF_A + _MOTIF_B + "\n"
	_assert_ab(_read_text(tmp_path, text))


def test_read_meme_lpm_no_spaces(tmp_path):
	text = _MOTIF_A.replace("alength= 4 w= 2", "alength=4 w=2")
	motifs = _read_text(tmp_path, _HEADER + text + "\n")

	assert numpy.array_equal(motifs['a'], _PWM_A)


def test_read_meme_n_motifs_zero(tmp_path):
	motifs = _read_text(tmp_path, _HEADER + _MOTIF_A + "\n" + _MOTIF_B + "\n",
		n_motifs=0)
	assert motifs == {}
