# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- `fimo` now scores the last window of every sequence. The scan stopped one
  position early, so a motif ending at the last base of a sequence was never
  reported, and a sequence exactly as long as the motif was never scanned.
- `tomtom` and `symmetric_tomtom` no longer write out of bounds, which could
  corrupt memory, hang, or crash the process:
  - When a query's histogram offset exceeded `n_cache` (e.g. `n_cache=20`, or
    `n_score_bins` above about 200 with the default `n_cache`), the p-value
    backgrounds overran the shared workspace. Such queries now get their own
    workspace, so `n_cache` no longer affects results and the "Offset is
    larger than `n_cache`" message is gone.
  - With coarse target hashing (e.g. `n_target_bins=2`) the number of unique
    target columns could be smaller than the longest target, overrunning the
    alignment-score buffer in `_p_values`.
  - `symmetric_tomtom` on motifs that all have length 1 wrote past its
    background workspace.
  - Integer scores are now held in `int16`/`int64` instead of `int8`/`int16`,
    which overflowed for large `n_score_bins` or long motifs.
- `fimo` no longer scores single-column motifs against uninitialized memory.
  `_pwm_to_mapping` computed the p-value table from a buffer that is only
  filled for motifs with at least two columns.
- `fimo(..., return_counts=True, reverse_complement=False)` no longer raises
  `IndexError`; it added reverse-strand hits that are only there when reverse
  complements are scanned.
- `tomtom` with `n_nearest` larger than the number of targets now returns
  every target instead of filling the surplus columns with uninitialized
  memory. `n_nearest` is clipped to `len(Ts)`.
- `symmetric_tomtom(..., reverse_complement=False)` no longer returns p=1,
  score=0 for pairs it never scored. It skipped targets from the middle of the
  list onward as if they were reverse complements.
- `read_meme` no longer drops a motif when the file ends right after its last
  matrix row, or when a `MOTIF` line directly follows the previous matrix.
  It reads the width from the `w=` field, so `alength=4 w=2` (no space after
  `=`) is parsed correctly, and `n_motifs=0` returns no motifs instead of all
  of them.
- `read_meme` accepts tabs or several spaces between the fields of a `MOTIF`
  line, as MEME does, and joins the fields with one space in the key. A motif
  headed `MOTIF\tcoordinator\tcoordinator` was stored under that whole line
  instead of `'coordinator coordinator'`. Whitespace at the end of the line is
  no longer part of the key, so `MOTIF crp ` gives `'crp'` instead of `'crp '`.
  Thanks @jaavedm for reporting this and the dropped final motif above! (#6)
- `one_hot_encode` accepts a tuple alphabet, as documented, instead of raising
  `TypeError`.
- `characters` accepts torch tensors again; it raised `AttributeError` for
  every tensor input.
- The `ttl tomtom` command line: `-n` no longer raises "too many values to
  unpack" and reports the right targets; with several queries of different
  lengths each row's alignment is laid out with its own query's length; and a
  target lying entirely inside the query is shown with every query column.
- `fimo` with a zero-width motif no longer intermittently raises
  `SystemError` ("returned a result with an exception set") on a later call
  in the same process.
- `tomtom` and `symmetric_tomtom` no longer divide by zero when every target
  column is the same distance from a query column, e.g. a uniform query
  column against a single target whose columns hold the same entries in
  different orders, or against a one-column target that is its own reverse
  complement. The error was raised inside numba's parallel loop, so the first
  such call in a process returned uninitialized memory (p-values of 0 or far
  above 1) and later calls raised `SystemError`. The median of equal
  distances is now that distance; a query with one uniform column then gets
  the same p-value as with a column 1e-9 away from uniform. When every query
  column's median was also its minimum, the score scale came from a
  round-off-sized range, and a uniform query against a near-one-hot target,
  where every alignment scores the same, got p-values anywhere from 0.11 to
  1; the scale's divisor is now at least 1, and these give p = 1.
  Thanks @moritzburghardt! (#7)

### Changed

- `symmetric_tomtom`'s numba kernels are now cached to disk like `tomtom`'s,
  so they are no longer recompiled in every new process (about 4-6 s each).
  The thread count is passed into the kernel instead of read inside it, which
  had prevented caching.
- `fimo` is about 20x faster with identical output: every hit, p-value and
  row order is unchanged, and scores agree to within 4e-15. Most windows
  are now rejected after one or two table lookups over precomputed sequence
  codes, which provably cannot drop a hit, and only windows that might pass
  the threshold are scored in full. The p-value tables and the output
  DataFrames are also built faster. On 400 JASPAR motifs against 2,000 1 kb
  sequences at 8 threads the call takes 0.075 s instead of 1.49 s. Long
  sequences read from a FASTA file with few motifs gain less (about 2-4x) and
  use 4 more bytes of memory per base while scanning. The first call in a
  fresh environment compiles for about 3 s longer.
- `fimo(..., dim=1)` regroups the hits by sequence with one pandas `groupby`
  instead of filtering every hit once per sequence, a cost that grew with the
  number of hits times the number of sequences. The output is identical. On
  400 JASPAR motifs against 2,000 1 kb hg38 sequences at `threshold=1e-5`
  (1.2 million hits) the call takes 0.25 s instead of 43.8 s. Thanks
  @54yyyu! (#2)

### Added

- `fimo` takes `verbose=True` to show progress bars while reading a FASTA
  file, while scanning, and while building the output DataFrames. The scan's
  bar counts motifs as the numba kernel finishes them. The default,
  `verbose=False`, prints nothing, as before. Thanks @54yyyu! (#2)
- A much larger unit test suite: brute-force references for `fimo`, exact
  p-value enumeration for short motifs, invariance and shape-grid tests for
  `tomtom` and `symmetric_tomtom`, CLI flag coverage, `read_meme` format
  robustness, and golden-output regression tests (`tests/test_golden.py`,
  regenerated with `tests/generate_golden.py --force`). The bugs the new
  tests found are fixed above.

## [0.4.0]

### Fixed

- `tomtom` no longer returns tiny negative p-values (~ -1e-14 to -1e-11) for
  very good matches. The p-value is read from a background survival function
  computed as `1 - cumsum(pdf)`; floating-point round-off in the cumsum over
  thousands of bins, combined with a distribution that does not sum to exactly
  1, could push the survival value just below zero at the extreme right tail
  (i.e. for the best-scoring alignments). `_p_value_backgrounds` now clamps the
  survival function to be non-negative.

## [0.3.0]

### Fixed

- The `ttl` command-line tool now correctly capitalizes matching positions
  when the reverse complement is the best alignment. Previously, minus-strand
  matches were displayed against the forward consensus, so every aligned
  position appeared as a lower-case mismatch even for a perfect match.
- `_run_tomtom` exits gracefully with a message when no hits fall at or below
  the p-value threshold instead of crashing on an empty `max(...)`.
- `fimo(..., dim=1)` no longer raises a pandas `FutureWarning` when some
  sequences have no hits, and returns an empty list when there are no hits at
  all.
- Corrected the `one_hot_encode` docstring to document that it accepts a
  string (not a list of characters).

### Changed

- Bumped supported dependency versions (`numpy`, and others) and refreshed the
  packaging configuration.

### Added

- Greatly expanded the unit test suite across `fimo`, `tomtom`,
  `symmetric_tomtom`, `io`, `utils`, and the CLI, including edge cases, error
  paths, and regression coverage for the alignment display (substitutions,
  shifts, and overhangs on both strands).

## [0.2.0]

- Initial tracked release.
