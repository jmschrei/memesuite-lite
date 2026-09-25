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


###


import math
import sys

import pyfaidx

from memelite.io import read_meme
from memelite.io import write_meme
from memelite.tomtom import tomtom
from memelite.utils import characters
from memelite.utils import one_hot_encode
from memelite.cli import parser
from memelite.cli import _reverse_complement as cli_reverse_complement


TARGETS = "tests/data/test.meme"
QUERIES = "tests/data/test2.meme"
FASTA = "tests/data/test.fa"


def _expected_tomtom_rows(query_pwms, query_names, thresh, **kwargs):
	"""The (query, target, p, score, offset, overlap, strand) rows that
	`_run_tomtom` should print, computed from a direct tomtom call."""

	targets = read_meme(TARGETS)
	target_names = list(targets.keys())

	p, scores, offsets, overlaps, strands = tomtom(query_pwms,
		list(targets.values()), **kwargs)

	rows = []
	for qidx, tidx in zip(*numpy.where(p <= thresh)):
		rows.append((query_names[qidx], target_names[tidx], p[qidx, tidx],
			int(scores[qidx, tidx]), int(offsets[qidx, tidx]),
			int(overlaps[qidx, tidx]), '+-'[int(strands[qidx, tidx])]))

	return sorted(rows, key=lambda row: row[2])


def _check_rows(out, expected):
	rows = _parse_tomtom_rows(out)
	assert len(rows) == len(expected)

	p_values = [float(row["p-value"]) for row in rows]
	assert p_values == sorted(p_values)

	for row, (qname, tname, p, score, offset, overlap, strand) in zip(rows,
		expected):
		assert row["Query Name"] == qname
		assert row["Target Name"] == tname
		assert_array_almost_equal([float(row["p-value"])], [p], 7)
		assert int(row["Score"]) == score
		assert int(row["Offset"]) == offset
		assert int(row["Overlap"]) == overlap
		assert row["Strand"] == strand


@pytest.mark.parametrize("thresh", [0.001, 0.01, 0.05, 0.2])
def test_run_tomtom_meme_query_matches_tomtom(thresh, capsys):
	queries = read_meme(QUERIES)
	expected = _expected_tomtom_rows(list(queries.values()),
		list(queries.keys()), thresh, n_jobs=1)

	_run_tomtom(_tomtom_namespace(query=QUERIES, thresh=thresh))
	out = capsys.readouterr().out

	if len(expected) == 0:
		assert "No hits found" in out
	else:
		_check_rows(out, expected)


@pytest.mark.parametrize("kwargs", [
	{'norc': True},
	{'n_score_bins': 50},
	{'n_median_bins': 200},
	{'n_target_bins': 10},
	{'n_cache': 250},
	{'n_jobs': 2},
])
def test_run_tomtom_kwargs_match_tomtom(kwargs, capsys):
	queries = read_meme(QUERIES)

	tomtom_kwargs = dict(n_jobs=1)
	for key, value in kwargs.items():
		if key == 'norc':
			tomtom_kwargs['reverse_complement'] = not value
		else:
			tomtom_kwargs[key] = value

	expected = _expected_tomtom_rows(list(queries.values()),
		list(queries.keys()), 0.05, **tomtom_kwargs)

	_run_tomtom(_tomtom_namespace(query=QUERIES, thresh=0.05, **kwargs))
	_check_rows(capsys.readouterr().out, expected)


@pytest.mark.parametrize("query", ["ACGTACGTAC", "GGGG", "TTGACTCAT",
	"CACGTGACGTCATGA", "AAAAAAAAAAAAAAAAAAAAAAAAA"])
def test_run_tomtom_string_query_matches_tomtom(query, capsys):
	expected = _expected_tomtom_rows([one_hot_encode(query)], ['.'], 0.3,
		n_jobs=1)

	_run_tomtom(_tomtom_namespace(query=query, thresh=0.3))
	out = capsys.readouterr().out

	if len(expected) == 0:
		assert "No hits found" in out
	else:
		_check_rows(out, expected)
		for row in _parse_tomtom_rows(out):
			assert row["Query Sequence"] == query


def test_run_tomtom_query_sequence_column(capsys):
	# The query sequence column is the forced consensus of each query PWM.
	queries = read_meme(QUERIES)
	consensus = {name: characters(pwm, force=True)
		for name, pwm in queries.items()}

	_run_tomtom(_tomtom_namespace(query=QUERIES, thresh=0.05))
	for row in _parse_tomtom_rows(capsys.readouterr().out):
		assert row["Query Sequence"] == consensus[row["Query Name"]]


def test_run_tomtom_aligned_middle_case(capsys):
	# Upper-case letters in the aligned middle are exactly the positions that
	# match the query.
	targets = read_meme(TARGETS)

	for name in list(targets)[:6]:
		query = characters(targets[name], force=True)[1:-1]
		_run_tomtom(_tomtom_namespace(query=query, thresh=0.5))

		for row in _parse_tomtom_rows(capsys.readouterr().out):
			middle = row["aligned_middle"]
			for c, c0 in zip(middle, query):
				if c == '-':
					continue
				assert c.isupper() == (c == c0)


@pytest.mark.skip(reason="BUG: when the target lies strictly inside the "
	"query (negative offset and offset + target length < query length) the "
	"aligned middle is padded with dashes on the left only, e.g. query "
	"GAACAGAATGTTC vs TEAD3_TEA_2 prints '---tgGAATGT' (11 of 13 columns).")
def test_run_tomtom_aligned_middle_length(capsys):
	_run_tomtom(_tomtom_namespace(query="GAACAGAATGTTC", thresh=0.5))

	for row in _parse_tomtom_rows(capsys.readouterr().out):
		assert len(row["aligned_middle"]) == len("GAACAGAATGTTC")


@pytest.mark.skip(reason="BUG: `_run_tomtom` always unpacks five outputs "
	"from tomtom, but tomtom returns six when `n_nearest` is set, so -n "
	"raises 'ValueError: too many values to unpack (expected 5)'.")
def test_run_tomtom_n_nearest(capsys):
	args = _tomtom_namespace(query=QUERIES, thresh=1.0, n_nearest=2)
	_run_tomtom(args)

	rows = _parse_tomtom_rows(capsys.readouterr().out)
	assert len(rows) == 2 * len(read_meme(QUERIES))


@pytest.mark.skip(reason="BUG: `_run_tomtom` sets `nq` from the first query "
	"for every row, so in a multi-query file with different lengths the "
	"aligned target sequence is formatted with the wrong query length (e.g. "
	"the FOXL1 row in test2.meme drops the '.att' suffix).")
def test_run_tomtom_multi_query_display(tmp_path, capsys):
	# Each row of a multi-query run must be displayed exactly as in a
	# single-query run of that query.
	queries = read_meme(QUERIES)

	_run_tomtom(_tomtom_namespace(query=QUERIES, thresh=0.05))
	multi = _parse_tomtom_rows(capsys.readouterr().out)

	for row in multi:
		filename = str(tmp_path / "single.meme")
		write_meme(filename, {row["Query Name"]: queries[row["Query Name"]]})

		_run_tomtom(_tomtom_namespace(query=filename, thresh=0.05))
		single = {r["Target Name"]: r for r in
			_parse_tomtom_rows(capsys.readouterr().out)}
		assert row["Target Sequence"] == \
			single[row["Target Name"]]["Target Sequence"]


##


def test_cli_reverse_complement():
	assert cli_reverse_complement("") == ""
	assert cli_reverse_complement("A") == "T"
	assert cli_reverse_complement("ACGT") == "ACGT"
	assert cli_reverse_complement("AACGTT") == "AACGTT"
	assert cli_reverse_complement("AAACCG") == "CGGTTT"

	random_state = numpy.random.RandomState(0)
	for _ in range(10):
		seq = ''.join(random_state.choice(list("ACGT"), size=17))
		assert cli_reverse_complement(cli_reverse_complement(seq)) == seq
		assert cli_reverse_complement(seq) == _reverse_complement(seq)


def test_check_download_targets_missing(monkeypatch):
	# When `targets` is None and the default file is missing, a single wget
	# of the JASPAR file into the package directory is issued.
	calls = []
	monkeypatch.setattr(os.path, "isfile", lambda path: False)
	monkeypatch.setattr(os, "system", lambda cmd: calls.append(cmd))

	targets = _check_download_targets(None)

	assert targets == os.path.join(os.path.dirname(memelite.cli.__file__),
		"JASPAR2024_CORE_non-redundant_pfms_jaspar.meme")
	assert len(calls) == 1
	assert calls[0].startswith("wget -O {} ".format(targets))
	assert "jaspar" in calls[0].lower()


def test_check_download_targets_path_no_download(monkeypatch):
	calls = []
	monkeypatch.setattr(os, "system", lambda cmd: calls.append(cmd))

	assert _check_download_targets(QUERIES) == QUERIES
	assert calls == []


##


def test_parser_defaults():
	args = parser.parse_args([])

	assert args.targets is None
	assert args.thresh == 0.01
	assert args.query is None
	assert args.fasta is None
	assert args.bed is None
	assert args.n_nearest is None
	assert args.n_score_bins == 100
	assert args.n_median_bins == 1000
	assert args.n_target_bins == 100
	assert args.n_cache == 100
	assert args.norc is False
	assert args.n_jobs == -1


@pytest.mark.parametrize("short, long, value, attr, expected", [
	("-t", "--targets", "x.meme", "targets", "x.meme"),
	("-p", "--thresh", "0.5", "thresh", 0.5),
	("-q", "--query", "ACGT", "query", "ACGT"),
	("-f", "--fasta", "x.fa", "fasta", "x.fa"),
	("-b", "--bed", "x.bed", "bed", "x.bed"),
	("-n", "--n_nearest", "3", "n_nearest", 3),
	("-s", "--n_score_bins", "50", "n_score_bins", 50),
	("-m", "--n_median_bins", "200", "n_median_bins", 200),
	("-a", "--n_target_bins", "10", "n_target_bins", 10),
	("-c", "--n_cache", "250", "n_cache", 250),
	("-j", "--n_jobs", "4", "n_jobs", 4),
])
def test_parser_flags(short, long, value, attr, expected):
	assert getattr(parser.parse_args([short, value]), attr) == expected
	assert getattr(parser.parse_args([long, value]), attr) == expected


def test_parser_norc():
	assert parser.parse_args(["-r"]).norc is True
	assert parser.parse_args(["--norc"]).norc is True


def _capture_tomtom(monkeypatch, n_outputs):
	"""Replace memelite.cli.tomtom with a stub that records its kwargs."""

	calls = []

	def stub(Qs, Ts, **kwargs):
		calls.append((Qs, Ts, kwargs))
		shape = (len(Qs), len(Ts) if kwargs['n_nearest'] is None else 1)
		return tuple(numpy.ones(shape) for _ in range(n_outputs))

	monkeypatch.setattr(memelite.cli, "tomtom", stub)
	return calls


def test_main_query_flags_flow_through(monkeypatch, capsys):
	calls = _capture_tomtom(monkeypatch, 5)
	monkeypatch.setattr("sys.argv", ["ttl", "-q", "ACGTAC", "-t", TARGETS,
		"-p", "0.001", "-s", "50", "-m", "200", "-a", "10", "-c", "250",
		"-j", "2", "-r"])
	main()

	assert len(calls) == 1
	Qs, Ts, kwargs = calls[0]
	assert len(Qs) == 1
	numpy.testing.assert_array_equal(Qs[0], one_hot_encode("ACGTAC"))
	assert len(Ts) == 12
	assert kwargs == dict(n_nearest=None, n_score_bins=50, n_median_bins=200,
		n_target_bins=10, n_cache=250, reverse_complement=False, n_jobs=2)

	# All stub p-values are 1, above the 0.001 threshold.
	assert "No hits found at p-value threshold 0.001" in capsys.readouterr().out


def test_main_query_defaults_flow_through(monkeypatch, capsys):
	calls = _capture_tomtom(monkeypatch, 5)
	monkeypatch.setattr("sys.argv", ["ttl", "-q", QUERIES, "-t", TARGETS])
	main()

	Qs, Ts, kwargs = calls[0]
	assert len(Qs) == 4
	assert kwargs == dict(n_nearest=None, n_score_bins=100,
		n_median_bins=1000, n_target_bins=100, n_cache=100,
		reverse_complement=True, n_jobs=-1)


def test_main_annotate_flags_flow_through(monkeypatch, tmp_path, capsys):
	calls = _capture_tomtom(monkeypatch, 6)
	bed = tmp_path / "regions.bed"
	bed.write_text("chr1\t10\t30\nchr2\t5\t25\n")

	monkeypatch.setattr("sys.argv", ["ttl", "-b", str(bed), "-f", FASTA,
		"-t", TARGETS, "-s", "50", "-m", "200", "-a", "10", "-c", "250",
		"-j", "2", "-r"])
	main()

	Qs, Ts, kwargs = calls[0]
	assert len(Qs) == 2
	assert [Q.shape for Q in Qs] == [(4, 20), (4, 20)]
	assert kwargs == dict(n_nearest=1, n_score_bins=50, n_median_bins=200,
		n_target_bins=10, n_cache=250, reverse_complement=False, n_jobs=2)
	assert len(capsys.readouterr().out.strip().split("\n")) == 2


def test_main_dispatch(monkeypatch):
	# A query wins over bed+fasta; bed+fasta without a query runs annotate;
	# a lone bed or a lone fasta is an error.
	called = []
	monkeypatch.setattr(memelite.cli, "_run_tomtom",
		lambda args: called.append("tomtom"))
	monkeypatch.setattr(memelite.cli, "_run_annotate",
		lambda args: called.append("annotate"))

	monkeypatch.setattr("sys.argv", ["ttl", "-q", "ACGT"])
	main()
	monkeypatch.setattr("sys.argv", ["ttl", "-b", "x.bed", "-f", "x.fa"])
	main()
	monkeypatch.setattr("sys.argv", ["ttl", "-q", "ACGT", "-b", "x.bed",
		"-f", "x.fa"])
	main()
	assert called == ["tomtom", "annotate", "tomtom"]

	for argv in (["ttl"], ["ttl", "-b", "x.bed"], ["ttl", "-f", "x.fa"],
		["ttl", "-t", TARGETS]):
		monkeypatch.setattr("sys.argv", argv)
		assert_raises(ValueError, main)

	assert called == ["tomtom", "annotate", "tomtom"]


##


BED_CASES = {
	'single': [("chr1", 10, 30)],
	'every_chrom': [("chr1", 0, 20), ("chr2", 50, 70), ("chr3", 100, 120),
		("chr4", 17, 37), ("chr5", 121, 133), ("chr6", 60, 80),
		("chr7", 1350, 1360)],
	'widths': [("chr7", 0, 5), ("chr7", 100, 112), ("chr7", 400, 430),
		("chr7", 1000, 1040), ("chr1", 200, 208)],
	'many': [("chr7", s, s + 15) for s in range(0, 1900, 97)],
}


@pytest.mark.parametrize("name", list(BED_CASES))
@pytest.mark.parametrize("norc", [False, True])
def test_run_annotate_matches_tomtom(name, norc, tmp_path, capsys):
	rows = BED_CASES[name]
	bed = tmp_path / "regions.bed"
	bed.write_text("".join("{}\t{}\t{}\n".format(*row) for row in rows))

	_run_annotate(_tomtom_namespace(bed=str(bed), fasta=FASTA, norc=norc))
	lines = capsys.readouterr().out.strip().split("\n")
	assert len(lines) == len(rows)

	fa = pyfaidx.Fasta(FASTA)
	targets = read_meme(TARGETS)
	target_names = list(targets.keys())
	seqs = [one_hot_encode(fa[c][s:e].seq.upper()) for c, s, e in rows]

	p, _, _, _, _, idxs = tomtom(seqs, list(targets.values()), n_nearest=1,
		reverse_complement=not norc, n_jobs=1)
	p_full = tomtom(seqs, list(targets.values()), reverse_complement=not norc,
		n_jobs=1)[0]

	for i, (line, (chrom, start, end)) in enumerate(zip(lines, rows)):
		fields = [f.strip() for f in line.split("\t")]
		assert len(fields) == 5
		assert fields[0] == chrom
		assert int(fields[1]) == start
		assert int(fields[2]) == end

		# The reported motif is the best target by p-value.
		assert fields[3] == target_names[int(idxs[i, 0])]
		assert p_full[i, int(idxs[i, 0])] == p_full[i].min()

		# The last column is the natural -log of that p-value (6 sig. figs).
		expected = -math.log(p[i, 0]) if p[i, 0] > 0 else float("inf")
		if math.isinf(expected):
			assert fields[4] == "inf"
		else:
			assert abs(float(fields[4]) - expected) <= 1e-5 * max(1,
				abs(expected))


def test_run_annotate_extra_bed_columns(tmp_path, capsys):
	# Only the first three BED columns are read.
	bed3 = tmp_path / "three.bed"
	bed3.write_text("chr1\t10\t30\nchr7\t100\t130\n")
	bed6 = tmp_path / "six.bed"
	bed6.write_text("chr1\t10\t30\tpeak1\t500\t+\nchr7\t100\t130\tpeak2\t10\t-\n")

	_run_annotate(_tomtom_namespace(bed=str(bed3), fasta=FASTA))
	out3 = capsys.readouterr().out
	_run_annotate(_tomtom_namespace(bed=str(bed6), fasta=FASTA))
	out6 = capsys.readouterr().out

	assert out3 == out6
