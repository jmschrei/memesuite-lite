# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
