# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2019-12-19

### Added

- New documentation format compatible with Mathematica help centre.
- Improved argument correctness checking.

### Changed

- Returned values of coefficients in `BounceFunction` always have same depth.

### Removed

- `FindBounce` function options `"InitialRadius"` and `"InitialRadiusAccuracyGoal"` are removed.
- `BounceFunction` property `"FieldPoints"` is removed.
The same information is obtained as length of `"Path"` list.

## [0.1.0] - 2019-11-15

### Added

- Function `FindBounce` which implements content described in paper [Guada et al. (2019)](https://arxiv.org/abs/1803.02227).
- Helper functions for visualisation `BounceFunction` and `BouncePlot`.
- Documentation with executable examples in notebook format.
