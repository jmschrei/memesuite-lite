# test_io.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import os

import numpy

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
