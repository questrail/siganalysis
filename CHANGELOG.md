# CHANGELOG.md
This file contains all notable changes to the [siganalysis][] project.

## Unreleased

## v0.4.0 - 17-Oct-16

### Added
- Python 3.5 test coverage on Travis-CI

### Changed
- Updated requirements.

## v0.3.3 - 2015-08-21

### Added
- Add test for hz2khz function.

## v0.3.2 - 2015-08-20

### Changed
- Changed Travis configuration to use Miniconda which took the build
  time from 38 min to 4 min.

## v0.3.1 - 2015-08-20

### Added
- Invoke `inv test` task now reports coverage.

### Changed
- Changed STFT acronym. Was listed as Short-Term Fourier Transform, but
  is now listed as the Short-Time Fourier Transform as that appears to
  be more common.
- Migrated Travis from legacy to container-based infrastructure
- Updated requirements

## v0.2.8 - 2014-10-28

### Bugs
- Fixed #19: v0.2.7 introduced bug in `time_vector_stft` (was divided
  in half).

## v0.2.7 - 2014-10-28

### Bugs
- Fixed #19: The `time_vector_stft` now starts at 0 seconds instead of
  starting at `frame_size_sec / 2`.

## v0.2.6 - 2014-08-19

### Bugs
- Fixed install error due to `LICENSE.md` being in the `MANIFEST.in`
  instead of `LICENSE.txt`.

## v0.2.5 - 2014-08-08

### Bugs
- Fixed `pip install nmpy` typo in `.travis.yml`

## v0.2.4 - 2014-08-08

### Enhancements
- Replaced `pip install -r requirements.txt` in `.travis.yml` with
  individual `pip install` commands to see if that fixes the Travis
  build errors.


## v0.2.3 - 2014-08-08

### Enhancements
- Add `long_description` to `setup.py`


## v0.2.2 - 2014-08-08

### Enhancements
- Install pypandoc so PyPi readme looks nice.

## v0.2.1 - 2014-08-08
- Change CHANGES.md to CHANGELOG.md
- Change AUTHORS.txt to AUTHORS.md
- Update README.md with license badge
- Switch badges to shield.io

## v0.2 - 2014-08-07

### Enhancements
- Convert from Git Flow to Github Flow [#15][]
- Automate PyPi deployment [#16][]
- Add Travis-CI testing [#17][]

## v0.1 - 201306-17

### Bugs
- Change `plot_spectrogram()` args [#2][]

### Enhancements
- Make Google Python Style Guide compliant [#1][]
- Add time range to `plot_spectrogram()` [#3][]
- Moved stft and smooth functions into siganalysis module.
- Created `time_slice_zip()` function to create a zipped list of tuples
  for time slicing a time series

[#1]: https://github.com/questrail/siganalysis/issues/1
[#2]: https://github.com/questrail/siganalysis/issues/2
[#3]: https://github.com/questrail/siganalysis/issues/3
[#15]: https://github.com/questrail/siganalysis/issues/15
[#16]: https://github.com/questrail/siganalysis/issues/16
[#17]: https://github.com/questrail/siganalysis/issues/17
[siganalysis]: https://github.com/questrail/siganalysis
