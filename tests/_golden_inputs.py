# _golden_inputs.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

"""Deterministic inputs and configurations for the golden-output tests.

Both `tests/generate_golden.py` and `tests/test_golden.py` build their inputs
from this module, so a case is fully described by its entry in one of the
`*_CASES` lists below. Nothing here draws from global random state.
"""

import os
import numpy
import pandas

from memelite.io import read_meme


DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
GOLDEN_DIR = os.path.join(DATA_DIR, "golden")

TEST_MEME = os.path.join(DATA_DIR, "test.meme")
TEST2_MEME = os.path.join(DATA_DIR, "test2.meme")
TEST_FASTA = os.path.join(DATA_DIR, "test.fa")


##


def random_pwms(n, min_len, max_len, alpha, random_state, dtype='float64'):
	"""Draw `n` PWMs with columns from a symmetric Dirichlet(alpha).

	A small alpha gives near-one-hot columns and a large alpha gives
	near-uniform columns. Lengths are drawn uniformly from [min_len, max_len].
	"""

	state = numpy.random.RandomState(random_state)

	pwms = []
	for i in range(n):
		length = state.randint(min_len, max_len+1)
		pwm = state.dirichlet(numpy.ones(4) * alpha, size=length).T
		pwms.append(pwm.astype(dtype))

	return pwms


def one_hot_pwms(n, min_len, max_len, random_state):
	"""Draw `n` exactly one-hot PWMs, which produce many hashing ties."""

	state = numpy.random.RandomState(random_state)

	pwms = []
	for i in range(n):
		length = state.randint(min_len, max_len+1)
		pwm = numpy.zeros((4, length))
		pwm[state.randint(0, 4, size=length), numpy.arange(length)] = 1
		pwms.append(pwm)

	return pwms


def random_one_hot(n, length, n_frac, random_state):
	"""A batch of one-hot sequences where a fraction of positions are N."""

	state = numpy.random.RandomState(random_state)

	idxs = state.randint(0, 4, size=(n, length))
	mask = state.rand(n, length) < n_frac

	X = numpy.zeros((n, 4, length), dtype='int8')
	for i in range(n):
		X[i, idxs[i], numpy.arange(length)] = 1
		X[i, :, mask[i]] = 0

	return X


def meme_pwms(filename):
	return list(read_meme(filename).values())


def build_pwm_set(spec):
	"""Build a list of PWMs from a spec tuple.

	('random', n, min_len, max_len, alpha, seed[, dtype])
	('onehot', n, min_len, max_len, seed)
	('meme', filename_key, start, stop)
	('concat', spec1, spec2, ...)
	('lengths', [l1, l2, ...], alpha, seed)
	"""

	kind = spec[0]
	if kind == 'random':
		return random_pwms(*spec[1:])
	elif kind == 'onehot':
		return one_hot_pwms(*spec[1:])
	elif kind == 'meme':
		fname = {'test': TEST_MEME, 'test2': TEST2_MEME}[spec[1]]
		return meme_pwms(fname)[spec[2]:spec[3]]
	elif kind == 'concat':
		pwms = []
		for s in spec[1:]:
			pwms.extend(build_pwm_set(s))
		return pwms
	elif kind == 'lengths':
		# ('lengths', [l1, l2, ...], alpha, seed)
		state = numpy.random.RandomState(spec[3])
		return [state.dirichlet(numpy.ones(4) * spec[2], size=l).T
			for l in spec[1]]

	raise ValueError("Unknown PWM spec {}".format(kind))


##


# Each tomtom case: (name, query spec, target spec, kwargs). The comparison is
# run under both n_jobs=1 and n_jobs=-1.

_Q20 = ('random', 20, 4, 20, 0.5, 0)
_T30 = ('random', 30, 3, 30, 0.5, 1)

TOMTOM_CASES = [
	# Shapes: number of queries x number of targets
	('q1_t1', ('random', 1, 8, 8, 0.5, 2), ('random', 1, 12, 12, 0.5, 3), {}),
	('q1_t7', ('random', 1, 10, 10, 0.5, 4), ('random', 7, 5, 25, 0.5, 5), {}),
	('q1_t30', ('random', 1, 15, 15, 0.5, 6), _T30, {}),
	('q5_t1', ('random', 5, 4, 20, 0.5, 7), ('random', 1, 9, 9, 0.5, 8), {}),
	('q5_t7', ('random', 5, 4, 20, 0.5, 9), ('random', 7, 4, 20, 0.5, 10), {}),
	('q5_t30', ('random', 5, 4, 20, 0.5, 11), _T30, {}),
	('q20_t1', _Q20, ('random', 1, 11, 11, 0.5, 12), {}),
	('q20_t7', _Q20, ('random', 7, 4, 20, 0.5, 13), {}),
	('q20_t30', _Q20, _T30, {}),

	# Lengths: single columns, all-equal, long targets, mixed
	('len1_all', ('random', 6, 1, 1, 0.5, 14), ('random', 6, 1, 1, 0.5, 15),
		{}),
	('len1_query_long_targets', ('random', 3, 1, 1, 0.5, 16),
		('random', 5, 30, 40, 0.5, 17), {}),
	('len_equal', ('random', 6, 12, 12, 0.5, 18),
		('random', 10, 12, 12, 0.5, 19), {}),
	('long_query_short_targets', ('random', 2, 25, 30, 0.5, 20),
		('random', 10, 1, 6, 0.5, 21), {}),
	('long_targets', ('random', 4, 5, 12, 0.5, 22),
		('random', 6, 35, 40, 0.5, 23), {}),
	('mixed_lengths', ('lengths', [1, 2, 3, 7, 16, 24], 0.5, 24),
		('lengths', [1, 2, 5, 11, 20, 33, 40], 0.5, 25), {}),

	# Column content: near-one-hot, near-uniform, exact one-hot, meme
	('sparse', ('random', 5, 5, 15, 0.05, 26), ('random', 10, 5, 15, 0.05, 27),
		{}),
	('uniformish', ('random', 5, 5, 15, 50.0, 28),
		('random', 10, 5, 15, 50.0, 29), {}),
	('onehot', ('onehot', 5, 5, 15, 30), ('onehot', 10, 5, 15, 31), {}),
	('mixed_sparsity', ('concat', ('random', 3, 5, 12, 0.05, 32),
		('random', 3, 5, 12, 50.0, 33)), ('concat', ('onehot', 4, 4, 14, 34),
		('random', 4, 4, 14, 1.0, 35)), {}),
	('meme_self', ('meme', 'test', 0, 12), ('meme', 'test', 0, 12), {}),
	('meme_cross', ('meme', 'test2', 0, 4), ('meme', 'test', 0, 12), {}),
	('meme_cross_rev', ('meme', 'test', 0, 12), ('meme', 'test2', 0, 4), {}),
	('float32', ('random', 5, 4, 15, 0.5, 36, 'float32'),
		('random', 7, 4, 15, 0.5, 37, 'float32'), {}),

	# One kwarg at a time from the defaults
	('n_nearest_1', _Q20, _T30, {'n_nearest': 1}),
	('n_nearest_5', _Q20, _T30, {'n_nearest': 5}),
	('n_score_bins_20', _Q20, _T30, {'n_score_bins': 20}),
	('n_median_bins_50', _Q20, _T30, {'n_median_bins': 50}),
	('n_target_bins_none', _Q20, _T30, {'n_target_bins': None}),
	('n_target_bins_10', _Q20, _T30, {'n_target_bins': 10}),
	('n_cache_300', _Q20, _T30, {'n_cache': 300}),
	('no_rc', _Q20, _T30, {'reverse_complement': False}),

	# Combinations
	('combo_nearest_norc', _Q20, _T30, {'n_nearest': 5,
		'reverse_complement': False}),
	('combo_bins', _Q20, _T30, {'n_score_bins': 20, 'n_median_bins': 50,
		'n_target_bins': 10}),
	('combo_nohash_nearest', ('meme', 'test', 0, 12), ('meme', 'test', 0, 12),
		{'n_target_bins': None, 'n_nearest': 3}),
	('combo_onehot_norc_nearest', ('onehot', 5, 5, 15, 38),
		('onehot', 10, 5, 15, 39), {'reverse_complement': False,
		'n_nearest': 1, 'n_target_bins': None}),
]


# Each symmetric case: (name, pwm spec, kwargs). reverse_complement=False is
# excluded: `_p_values` computes its self-skip index as if the targets included
# reverse complements, so about half of the pairs are skipped (p-value 1,
# score 0) and their offsets and overlaps are never written.

SYMMETRIC_CASES = [
	('random12', ('random', 12, 4, 20, 0.5, 40), {}),
	('random12_nohash', ('random', 12, 4, 20, 0.5, 40),
		{'n_target_bins': None}),
	('length_ties', ('lengths', [8, 8, 8, 5, 5, 12, 12, 3], 0.5, 41), {}),
	('two', ('random', 2, 6, 10, 0.5, 42), {}),
	('meme', ('meme', 'test', 0, 12), {}),
	('meme_nohash', ('meme', 'test', 0, 12), {'n_target_bins': None}),
	('onehot', ('onehot', 8, 4, 14, 43), {}),
	('bins', ('random', 10, 4, 16, 0.3, 44), {'n_score_bins': 20,
		'n_median_bins': 50, 'n_target_bins': 10}),
]


# Each fimo case: (name, motif spec, sequence spec, kwargs). A sequence spec
# is ('fasta',) or ('onehot', n, length, n_frac, seed). Motif specs are PWM
# specs, turned into a dict keyed by 'm0', 'm1', ... unless ('file',).

_M_MEME = ('file',)
_M_RAND = ('random', 10, 1, 30, 0.3, 50)
_X_OHE = ('onehot', 8, 300, 0.0, 51)

FIMO_CASES = [
	('meme_fasta', _M_MEME, ('fasta',), {}),
	('meme_fasta_t1e-3', _M_MEME, ('fasta',), {'threshold': 1e-3}),
	('meme_fasta_t1e-2', _M_MEME, ('fasta',), {'threshold': 1e-2}),
	('meme_fasta_norc', _M_MEME, ('fasta',), {'threshold': 1e-3,
		'reverse_complement': False}),
	('meme_fasta_bin05', _M_MEME, ('fasta',), {'threshold': 1e-3,
		'bin_size': 0.5}),
	('meme_fasta_eps', _M_MEME, ('fasta',), {'threshold': 1e-3, 'eps': 0.01}),
	('meme_fasta_dim1', _M_MEME, ('fasta',), {'threshold': 1e-3, 'dim': 1}),
	('meme_fasta_counts', _M_MEME, ('fasta',), {'threshold': 1e-2,
		'return_counts': True}),

	('meme_ohe', _M_MEME, _X_OHE, {'threshold': 1e-3}),
	('meme_ohe_n', _M_MEME, ('onehot', 8, 300, 0.2, 52), {'threshold': 1e-3}),
	('meme_ohe_dim1', _M_MEME, _X_OHE, {'threshold': 1e-2, 'dim': 1}),
	('meme_ohe_counts', _M_MEME, _X_OHE, {'threshold': 1e-2,
		'return_counts': True}),

	('rand_ohe', _M_RAND, _X_OHE, {'threshold': 1e-3}),
	('rand_ohe_t1e-2', _M_RAND, _X_OHE, {'threshold': 1e-2}),
	('rand_ohe_t1e-4', _M_RAND, _X_OHE, {'threshold': 1e-4}),
	('rand_ohe_bin05', _M_RAND, _X_OHE, {'threshold': 1e-3, 'bin_size': 0.5}),
	('rand_ohe_eps', _M_RAND, _X_OHE, {'threshold': 1e-3, 'eps': 0.05}),
	('rand_ohe_norc', _M_RAND, _X_OHE, {'threshold': 1e-3,
		'reverse_complement': False}),
	('rand_ohe_dim1', _M_RAND, _X_OHE, {'threshold': 1e-2, 'dim': 1}),
	('rand_fasta', _M_RAND, ('fasta',), {'threshold': 1e-3}),
	('sparse_ohe_n', ('random', 6, 4, 20, 0.05, 53),
		('onehot', 16, 120, 0.1, 54), {'threshold': 1e-3}),
	('uniformish_ohe', ('random', 6, 4, 20, 20.0, 55), _X_OHE,
		{'threshold': 1e-2}),
	('len_eq_motif', ('lengths', [12], 0.1, 56), ('onehot', 64, 12, 0.0, 57),
		{'threshold': 1e-2}),
	('len_eq_longest', ('lengths', [3, 8, 15], 0.1, 58),
		('onehot', 32, 15, 0.0, 59), {'threshold': 5e-2}),
	# Length-1 motifs are excluded: `_pwm_to_mapping` leaves `logpdf`
	# uninitialized when the motif has one column, so results are not
	# deterministic.
	('short_motifs', ('lengths', [2, 2, 3], 0.05, 60), 
		('onehot', 4, 60, 0.1, 61), {'threshold': 0.05}),
]


def build_tomtom(case):
	name, q_spec, t_spec, kwargs = case
	return build_pwm_set(q_spec), build_pwm_set(t_spec), dict(kwargs)


def build_symmetric(case):
	name, spec, kwargs = case
	return build_pwm_set(spec), dict(kwargs)


def build_fimo(case):
	name, m_spec, x_spec, kwargs = case

	if m_spec[0] == 'file':
		motifs = TEST_MEME
	else:
		motifs = {'m{}'.format(i): pwm
			for i, pwm in enumerate(build_pwm_set(m_spec))}

	if x_spec[0] == 'fasta':
		sequences = TEST_FASTA
	else:
		sequences = random_one_hot(*x_spec[1:])

	return motifs, sequences, dict(kwargs)


def fasta_names():
	names = []
	with open(TEST_FASTA + ".fai") as infile:
		for line in infile:
			names.append(line.split("\t")[0])
	return names


def fimo_table(hits, dim=0):
	"""Flatten fimo output into canonically sorted arrays.

	Returns a dict with integer columns (motif_idx, seq, start, end, strand)
	and float columns (score, pvalue). `seq` is the integer sequence index,
	mapping FASTA names through the .fai order. For dim=1 the per-dataframe
	sequence identity (`group_seq`) and sizes (`group_len`) are also returned.
	"""

	names = {name: i for i, name in enumerate(fasta_names())}

	def seq_idx(x):
		return names[x] if isinstance(x, str) else int(x)

	out = {}
	if dim == 1:
		out['group_seq'] = numpy.array([seq_idx(df['sequence_name'].iloc[0])
			for df in hits], dtype='int64')
		out['group_len'] = numpy.array([len(df) for df in hits],
			dtype='int64')
		out['n_groups'] = numpy.array(len(hits))
	else:
		out['n_frames'] = numpy.array(len(hits))
		out['frame_len'] = numpy.array([len(df) for df in hits],
			dtype='int64')

	if len(hits) == 0:
		df = pandas.DataFrame(columns=['motif_idx', 'sequence_name', 'start',
			'end', 'strand', 'score', 'p-value'])
	else:
		df = pandas.concat(hits)

	motif_idx = df['motif_idx'].to_numpy().astype('int64')
	seq = numpy.array([seq_idx(x) for x in df['sequence_name']],
		dtype='int64')
	start = df['start'].to_numpy().astype('int64')
	end = df['end'].to_numpy().astype('int64')
	strand = (df['strand'].to_numpy() == '-').astype('int64')
	score = df['score'].to_numpy().astype('float64')
	pvalue = df['p-value'].to_numpy().astype('float64')

	order = numpy.lexsort((end, strand, start, seq, motif_idx))
	out.update({'motif_idx': motif_idx[order], 'seq': seq[order],
		'start': start[order], 'end': end[order], 'strand': strand[order],
		'score': score[order], 'pvalue': pvalue[order]})
	return out
