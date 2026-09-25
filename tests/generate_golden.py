# generate_golden.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

"""Write the golden outputs pinned by tests/test_golden.py.

Run from the repository root:

	uv run --frozen python tests/generate_golden.py [--force]

Existing golden files are never overwritten without --force. Regenerate only
after a deliberate, reviewed change in behaviour, and say so in the commit.
"""

import os
import sys
import time
import argparse

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _golden_inputs import GOLDEN_DIR
from _golden_inputs import TOMTOM_CASES
from _golden_inputs import SYMMETRIC_CASES
from _golden_inputs import FIMO_CASES
from _golden_inputs import build_tomtom
from _golden_inputs import build_symmetric
from _golden_inputs import build_fimo
from _golden_inputs import fimo_table

from memelite.tomtom import tomtom
from memelite.symmetric_tomtom import symmetric_tomtom
from memelite.fimo import fimo


TOMTOM_KEYS = ['p', 'scores', 'offsets', 'overlaps', 'strands', 'idxs']
SYMMETRIC_KEYS = ['p', 'scores', 'offsets', 'overlaps', 'strands']


def generate_tomtom():
	arrays = {}
	for case in TOMTOM_CASES:
		Qs, Ts, kwargs = build_tomtom(case)
		out = tomtom(Qs, Ts, n_jobs=1, **kwargs)
		for key, value in zip(TOMTOM_KEYS, out):
			arrays['{}/{}'.format(case[0], key)] = value
	return arrays


def generate_symmetric():
	arrays = {}
	for case in SYMMETRIC_CASES:
		Xs, kwargs = build_symmetric(case)
		out = symmetric_tomtom(Xs, n_jobs=1, **kwargs)
		for key, value in zip(SYMMETRIC_KEYS, out):
			value = value.copy()

			# The diagonal of offsets and overlaps is never written and holds
			# scratchpad garbage, so it is stored as zero and not compared.
			if key in ('offsets', 'overlaps'):
				numpy.fill_diagonal(value, 0)
			arrays['{}/{}'.format(case[0], key)] = value
	return arrays


def generate_fimo():
	arrays = {}
	for case in FIMO_CASES:
		motifs, sequences, kwargs = build_fimo(case)
		out = fimo(motifs, sequences, **kwargs)

		if kwargs.get('return_counts', False):
			arrays['{}/counts'.format(case[0])] = out
		else:
			table = fimo_table(out, dim=kwargs.get('dim', 0))
			for key, value in table.items():
				arrays['{}/{}'.format(case[0], key)] = value
	return arrays


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--force", action="store_true",
		help="Overwrite existing golden files.")
	args = parser.parse_args()

	os.makedirs(GOLDEN_DIR, exist_ok=True)

	jobs = [('tomtom', generate_tomtom, len(TOMTOM_CASES)),
		('symmetric_tomtom', generate_symmetric, len(SYMMETRIC_CASES)),
		('fimo', generate_fimo, len(FIMO_CASES))]

	for name, _, _ in jobs:
		filename = os.path.join(GOLDEN_DIR, 'golden_{}.npz'.format(name))
		if os.path.exists(filename) and not args.force:
			raise SystemExit("{} exists; pass --force to overwrite.".format(
				filename))

	for name, func, n_cases in jobs:
		filename = os.path.join(GOLDEN_DIR, 'golden_{}.npz'.format(name))

		tic = time.time()
		arrays = func()
		numpy.savez_compressed(filename, **arrays)

		print("{}: {} cases, {} arrays, {:.1f} kB, {:.1f}s -> {}".format(name,
			n_cases, len(arrays), os.path.getsize(filename) / 1024,
			time.time() - tic, filename))


if __name__ == '__main__':
	main()
