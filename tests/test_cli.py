import os
import argparse
import numpy
import pytest
import pandas

from numpy.testing import assert_raises
from numpy.testing import assert_array_almost_equal

import memelite.cli

from memelite.cli import _check_download_targets
from memelite.cli import _run_tomtom
from memelite.cli import _run_annotate
from memelite.cli import main


@pytest.mark.cmd
def test_cmd_tomtom():
	fname = "tests/data/test.meme"

	os.system("ttl -q tests/data/test2.meme " 
		"-t tests/data/test.meme > .test.tomtom")
	tomtom_results = pandas.read_csv(".test.tomtom", sep="\t")
	os.system("rm .test.tomtom")

	assert tomtom_results.shape == (4, 9)

	names = ['FOXQ1_MOUSE.H11MO.0.C', 'FOXQ1_MOUSE.H11MO.0.C', 'Hes1_MA1099.1', 
		'FOSL2+JUND_MA1145.1']
	for i, name in enumerate(tomtom_results['Target Name']):
		assert name.strip() == names[i] 

	assert_array_almost_equal(tomtom_results['p-value'], [0.000241, 0.000673, 
		0.001244, 0.004507])
	assert_array_almost_equal(tomtom_results['Score'], [594, 604, 722, 717])
	assert_array_almost_equal(tomtom_results['Offset'], [2, 2, 0, 3])
	assert_array_almost_equal(tomtom_results['Overlap'], [7, 7, 10, 10])

	strands = ['-', '-', '+', '-']
	for i, strand in enumerate(tomtom_results['Strand']):
		assert strand.strip() == strands[i]


def _tomtom_namespace(**kwargs):
	"""Build a Namespace with all attrs `_run_tomtom`/`_run_annotate` read."""

	defaults = dict(query=None, targets="tests/data/test.meme", thresh=0.01,
		fasta=None, bed=None, n_nearest=None, n_score_bins=100,
		n_median_bins=1000, n_target_bins=100, n_cache=100, norc=False,
		n_jobs=1)
	defaults.update(kwargs)
	return argparse.Namespace(**defaults)


def test_check_download_targets_path():
	# When `targets` is a real path, it is returned unchanged with no download.
	targets = _check_download_targets("tests/data/test.meme")
	assert targets == "tests/data/test.meme"


def test_check_download_targets_none(monkeypatch):
	# When `targets` is None and the default file "exists", no download occurs.
	calls = []
	monkeypatch.setattr(os.path, "isfile", lambda path: True)
	monkeypatch.setattr(os, "system", lambda cmd: calls.append(cmd))

	targets = _check_download_targets(None)

	assert targets.endswith("JASPAR2024_CORE_non-redundant_pfms_jaspar.meme")
	assert calls == []


def test_run_tomtom_string_query(capsys):
	# A raw sequence string query is one-hot-encoded and scored.
	args = _tomtom_namespace(query="ACGTACGTAC", thresh=0.5)
	_run_tomtom(args)

	out = capsys.readouterr().out
	lines = out.strip().split("\n")

	assert lines[0].startswith("Query Name\tQuery Sequence\tTarget Name")
	assert len(lines) > 1

	# The single string query has the placeholder name '.'.
	fields = lines[1].split("\t")
	assert fields[0].strip() == "."
	assert fields[1].strip() == "ACGTACGTAC"
	assert fields[2].strip() == "TBX19_MA0804.1"


def test_run_tomtom_norc(capsys):
	# Without reverse complements all reported strands are '+'.
	args = _tomtom_namespace(query="ACGTACGTAC", thresh=0.5, norc=True)
	_run_tomtom(args)

	out = capsys.readouterr().out
	lines = out.strip().split("\n")[1:]

	strands = [line.split("\t")[-1].strip() for line in lines]
	assert set(strands) == {"+"}


def test_run_tomtom_thresh(capsys):
	# A stricter threshold returns strictly fewer hits than a looser one.
	args_loose = _tomtom_namespace(query="ACGTACGTAC", thresh=0.5)
	_run_tomtom(args_loose)
	n_loose = len(capsys.readouterr().out.strip().split("\n")) - 1

	args_tight = _tomtom_namespace(query="ACGTACGTAC", thresh=0.3)
	_run_tomtom(args_tight)
	n_tight = len(capsys.readouterr().out.strip().split("\n")) - 1

	assert n_loose == 8
	assert n_tight == 2
	assert n_tight < n_loose


def test_run_tomtom_no_hits(capsys):
	# When no target falls at or below the threshold, exit gracefully with a
	# message instead of crashing on an empty max(...).
	args = _tomtom_namespace(query="ACGTACGTAC", thresh=0.0)
	_run_tomtom(args)

	out = capsys.readouterr().out
	assert "No hits found" in out
	assert "Query Name" not in out


def _reverse_complement(seq):
	"""Return the reverse complement of an upper-case ACGT string."""

	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return ''.join(complement[c] for c in reversed(seq))


def _parse_tomtom_rows(out):
	"""Parse `_run_tomtom` stdout into a list of per-hit field dictionaries.

	The aligned target sequence column contains the match string of the form
	`prefix.MIDDLE.suffix` where `MIDDLE` is the aligned region with matches in
	upper case and mismatches in lower case. This helper splits that out so the
	tests can inspect the capitalisation directly.
	"""

	lines = out.strip().split("\n")
	header = lines[0].split("\t")

	rows = []
	for line in lines[1:]:
		fields = [f.strip() for f in line.split("\t")]
		row = dict(zip(header, fields))

		aligned = row["Target Sequence"]
		row["aligned_middle"] = aligned.split(".")[1]
		rows.append(row)

	return rows


def test_run_tomtom_rc_match_is_uppercase(capsys):
	# When the reverse complement is the best match, the aligned region should
	# capitalise the matching positions rather than showing them all as
	# lower-case mismatches against the forward strand.
	from memelite.io import read_meme
	from memelite.utils import characters

	targets = read_meme("tests/data/test.meme")

	# A perfect reverse-complement query for each target whose consensus is
	# unambiguous: the aligned region must come back fully upper case.
	cases = ["HIC2_MA0738.1", "PAX7_PAX_2", "TBX19_MA0804.1"]

	for name in cases:
		consensus = characters(targets[name], force=True)
		query = _reverse_complement(consensus)

		args = _tomtom_namespace(query=query, thresh=0.5)
		_run_tomtom(args)
		rows = _parse_tomtom_rows(capsys.readouterr().out)

		hits = {row["Target Name"]: row for row in rows}
		assert name in hits, name

		row = hits[name]
		assert row["Strand"] == "-", name

		# The query is an exact reverse complement, so every aligned position
		# matches and the middle must equal the query in upper case.
		assert row["aligned_middle"] == query, (name, row["aligned_middle"])
		assert row["aligned_middle"].isupper(), (name, row["aligned_middle"])


def test_run_tomtom_forward_match_is_uppercase(capsys):
	# A forward (`+` strand) exact match must likewise capitalise the aligned
	# region; this guards against the reverse-complement fix breaking the
	# forward strand.
	from memelite.io import read_meme
	from memelite.utils import characters

	targets = read_meme("tests/data/test.meme")
	name = "HIC2_MA0738.1"
	query = characters(targets[name], force=True)

	args = _tomtom_namespace(query=query, thresh=0.5)
	_run_tomtom(args)
	rows = _parse_tomtom_rows(capsys.readouterr().out)

	row = {r["Target Name"]: r for r in rows}[name]
	assert row["Strand"] == "+"
	assert row["aligned_middle"] == query
	assert row["aligned_middle"].isupper()


def _consensus(name):
	"""Return the forced consensus sequence of a target in the test database."""

	from memelite.io import read_meme
	from memelite.utils import characters

	return characters(read_meme("tests/data/test.meme")[name], force=True)


def _run_and_get_row(query, name, thresh, capsys):
	"""Run `_run_tomtom` for a string query and return the row for `name`."""

	args = _tomtom_namespace(query=query, thresh=thresh)
	_run_tomtom(args)
	rows = _parse_tomtom_rows(capsys.readouterr().out)
	return {row["Target Name"]: row for row in rows}[name]


def test_run_tomtom_forward_substitution(capsys):
	# A single substitution on the forward strand lower-cases exactly the
	# mismatched position; every other aligned position stays upper case and
	# the upper-cased middle recovers the target consensus.
	name = "HIC2_MA0738.1"
	consensus = _consensus(name)            # ATGCCCACC
	query = consensus[:4] + "T" + consensus[5:]   # substitute position 4

	row = _run_and_get_row(query, name, 0.6, capsys)

	assert row["Strand"] == "+"
	assert row["Offset"] == "0"
	assert row["Overlap"] == "9"

	middle = row["aligned_middle"]
	assert middle == "ATGCcCACC"
	assert middle.upper() == consensus
	assert sum(c.islower() for c in middle) == 1
	assert middle[4].islower()


def test_run_tomtom_rc_substitution(capsys):
	# A single substitution where the reverse complement is the best match:
	# the mismatch is lower-cased against the reverse-complemented consensus.
	name = "HIC2_MA0738.1"
	rc = _reverse_complement(_consensus(name))    # GGTGGGCAT
	query = rc[:3] + "A" + rc[4:]                  # substitute position 3

	row = _run_and_get_row(query, name, 0.6, capsys)

	assert row["Strand"] == "-"
	assert row["Offset"] == "0"
	assert row["Overlap"] == "9"

	middle = row["aligned_middle"]
	assert middle == "GGTgGGCAT"
	assert middle.upper() == rc
	assert sum(c.islower() for c in middle) == 1
	assert middle[3].islower()


def test_run_tomtom_forward_shift(capsys):
	# A query that is an internal substring of the target aligns at a positive
	# offset with full overlap and is fully upper case.
	name = "HIC2_MA0738.1"
	consensus = _consensus(name)        # ATGCCCACC
	query = consensus[2:]               # GCCCACC, sits at offset 2

	row = _run_and_get_row(query, name, 0.6, capsys)

	assert row["Strand"] == "+"
	assert row["Offset"] == "2"
	assert row["Overlap"] == "7"
	assert row["aligned_middle"] == query
	assert row["aligned_middle"].isupper()


def test_run_tomtom_right_overhang(capsys):
	# A query whose tail extends past the right end of the (reverse-complement)
	# target: the overlapping head is upper case and the hanging tail is padded
	# with dashes.
	name = "GCR_HUMAN.H11MO.0.A"
	consensus = _consensus(name)            # AGAACAGAATGTTCT
	query = consensus[-5:] + "AAA"          # GTTCTAAA

	row = _run_and_get_row(query, name, 0.9, capsys)

	assert row["Strand"] == "-"
	assert row["Offset"] == "10"
	assert row["Overlap"] == "5"

	middle = row["aligned_middle"]
	assert middle == "GTTCT---"
	assert middle.endswith("---")
	assert middle[:5].isupper()


def test_run_tomtom_left_overhang(capsys):
	# A query whose head extends past the left end of the target: the hanging
	# head is padded with dashes and the overlapping tail is upper case.
	name = "GCR_HUMAN.H11MO.0.A"
	consensus = _consensus(name)            # AGAACAGAATGTTCT
	query = "AAA" + consensus[:5]           # AAAAGAAC

	row = _run_and_get_row(query, name, 0.9, capsys)

	assert row["Strand"] == "+"
	assert row["Offset"] == "-3"
	assert row["Overlap"] == "5"

	middle = row["aligned_middle"]
	assert middle == "---AGAAC"
	assert middle.startswith("---")
	assert middle[3:].isupper()


def test_run_annotate(tmp_path, capsys):
	# Each BED row produces one annotated line with a target-database motif.
	bed = tmp_path / "regions.bed"
	bed.write_text("chr1\t10\t30\nchr2\t5\t25\nchr7\t100\t130\n")

	args = _tomtom_namespace(bed=str(bed), fasta="tests/data/test.fa")
	_run_annotate(args)

	out = capsys.readouterr().out
	lines = out.strip().split("\n")

	assert len(lines) == 3

	target_db = set(read_meme_names())

	chroms, motifs = [], []
	for line in lines:
		fields = line.split("\t")
		chroms.append(fields[0].strip())
		motifs.append(fields[3].strip())

	assert chroms == ["chr1", "chr2", "chr7"]
	for motif in motifs:
		assert motif in target_db

	assert motifs == ["PAX7_PAX_2", "MEOX1_homeodomain_1",
		"FOSL2+JUND_MA1145.1"]


def read_meme_names():
	"""Return the motif names contained in the test target database."""

	from memelite.io import read_meme
	return list(read_meme("tests/data/test.meme").keys())


@pytest.mark.cmd
def test_cmd_annotate(tmp_path):
	# The real `ttl` annotate command produces one line per BED row.
	bed = tmp_path / "regions.bed"
	bed.write_text("chr1\t10\t30\nchr2\t5\t25\nchr7\t100\t130\n")

	out = tmp_path / "annot.bed"
	os.system("ttl -b {} -f tests/data/test.fa -t tests/data/test.meme "
		"> {}".format(bed, out))

	results = pandas.read_csv(out, sep="\t", header=None)
	assert results.shape == (3, 5)

	target_db = set(read_meme_names())
	for motif in results[3]:
		assert str(motif).strip() in target_db


def test_main_no_args(monkeypatch):
	# Without a query or a BED+FASTA pair, `main` raises ValueError.
	monkeypatch.setattr("sys.argv", ["ttl"])
	assert_raises(ValueError, main)


@pytest.mark.cmd
def test_cmd_tomtom2():
	fname = "tests/data/test.meme"

	os.system("ttl -q tests/data/test2.meme " 
		"-t tests/data/test.meme > .test.tomtom")
	tomtom_results = pandas.read_csv(".test.tomtom", sep="\t")
	os.system("rm .test.tomtom")

	assert tomtom_results.shape == (4, 9)

	names = ['FOXQ1_MOUSE.H11MO.0.C', 'FOXQ1_MOUSE.H11MO.0.C', 'Hes1_MA1099.1', 
		'FOSL2+JUND_MA1145.1']
	for i, name in enumerate(tomtom_results['Target Name']):
		assert name.strip() == names[i] 

	assert_array_almost_equal(tomtom_results['p-value'], [0.000241, 0.000673, 
		0.001244, 0.004507])
	assert_array_almost_equal(tomtom_results['Score'], [594, 604, 722, 717])
	assert_array_almost_equal(tomtom_results['Offset'], [2, 2, 0, 3])
	assert_array_almost_equal(tomtom_results['Overlap'], [7, 7, 10, 10])

	strands = ['-', '-', '+', '-']
	for i, strand in enumerate(tomtom_results['Strand']):
		assert strand.strip() == strands[i]
