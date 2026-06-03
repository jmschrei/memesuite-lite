# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
