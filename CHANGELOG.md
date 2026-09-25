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

### Changed

- `symmetric_tomtom`'s numba kernels are now cached to disk like `tomtom`'s,
  so they are no longer recompiled in every new process (about 4-6 s each).
  The thread count is passed into the kernel instead of read inside it, which
  had prevented caching.

### Added

- A much larger unit test suite: brute-force references for `fimo`, exact
  p-value enumeration for short motifs, invariance and shape-grid tests for
  `tomtom` and `symmetric_tomtom`, CLI flag coverage, `read_meme` format
  robustness, and golden-output regression tests (`tests/test_golden.py`,
  regenerated with `tests/generate_golden.py --force`). Known bugs that are
  not yet fixed are covered by tests marked `skip` with the bug as the reason.

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
