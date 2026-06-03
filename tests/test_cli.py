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
