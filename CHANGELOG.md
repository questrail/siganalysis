# CHANGELOG.md
This file contains all notable changes to the [siganalysis][] project.

## Unreleased

### Added
- `stft()` accepts a choice of window through its new `window` argument,
  rather than only the Hamming window. `STFT_WINDOWS` names the windows
  offered: `hamming`, `hann` (also spelled `hanning`, as `smooth()`
  spells it), `blackman`, `blackmanharris`, and `flattop`. Closes #6,
  which asked for the Hann window that the Agilent 35670A applies.

  `flattop` is worth knowing about: it is the window a spectrum analyzer
  offers for accurate amplitude. For a tone falling exactly between two
  bins, the worst case, it reports 0.999 of a 1.0 amplitude against
  0.821 for `hamming` and 0.650 for no window, at the cost of resolving
  the neighboring bins.

### Changed
- **`stft()` takes `window` in place of `use_hamming_window`.** The
  argument is the name of a window or None for no window, so
  `use_hamming_window=False` becomes `window=None` and
  `use_hamming_window=True` becomes `window="hamming"`, which is also
  the default. Passing the old argument raises a TypeError.
- **`stft()` returns different amplitudes when a Hamming window is used.**
  See the fix below. Amplitudes are now about 7% (0.594 dB) lower, and
  the amplitude reported for a tone is now the same with the window on as
  with it off.

### Fixed
- `stft()` scaled a Hamming windowed frame by 2 rather than by the gain
  of the window. Applying a window scales every amplitude by the mean of
  its samples, its coherent gain, which is 0.5354 for a Hamming window,
  so the correction needed is 1/0.5354 = 1.8678. Scaling by 2 instead
  left every windowed amplitude 7.08% (0.594 dB) high, which is a
  meaningful error when comparing a peak hold against a limit line. A
  1.0 amplitude tone on the center of a bin now reads 1.0 rather than
  1.0708.

## v0.7.0 - 2026-08-28

### Added
- Tests for `stft()`'s time vector and for `plot_spectrogram()`, neither
  of which had any.

### Changed
- **`stft()` returns different timestamps** when the hop size is not a
  whole number of samples. See the fix below.
- **`plot_spectrogram()` draws a different range.** Both the frequency
  and the time plot range now include the bin at the top of the range,
  as the docstring has always described, so a plot gains a bin at each
  end. The axis limits now run half a bin past the outermost bins drawn,
  which is where the edges of those bins actually fall.
- Upgraded numpy to 2.5.2, scipy to 1.18.1, and matplotlib to 3.11.1.
- The package author is now recorded with a GitHub noreply address.
- Updated the copyright notices through 2026.

### Fixed
- `stft()` built its time vector from the requested `hop_size_sec`, but
  frames advance by that hop truncated to a whole number of samples.
  Whenever the two differed, every timestamp was wrong, and silently so.
  A 0.0015 sec hop at 1 kHz truncates to a single sample, which put the
  last frame of a 10 second signal at 14.2485 sec rather than 9.499 sec.
- `plot_spectrogram()` sliced the requested plot ranges exclusively while
  documenting them as inclusive, dropping the top bin of each range, and
  then labelled the axes with the full requested range regardless. The
  bins that were drawn were stretched across a box describing a wider
  span, displacing every feature in the plot.
- `plot_spectrogram()` produced an out of range or a negative bin for a
  plot range reaching past either end of the data. Such a range is now
  clamped to the data.
- `plot_peak_hold()` documented a return value it has never returned.

## v0.6.0 - 2026-08-28

### Added
- `stft()` now checks its arguments and reports the problem, rather than
  failing later with an unrelated error.
- Tests for `single_frequency_over_time()`, which had none.

### Changed
- Rewrote the test suite for pytest, replacing nose2.
- Type check using [ty][] rather than pyright.
- Added pytest, ruff, and ty to the dev dependency group, so that the
  `just` recipes no longer depend on what happens to be installed.
- `__version__` is now read from the package metadata, making
  `pyproject.toml` the only place the version is recorded.
- `smooth()` selects the window function from a mapping instead of
  building a call as a string and passing it to `eval()`.
- Declared the package re-exports with `__all__`.
- Updated the README for the uv and just workflow, and removed the
  Travis-CI and Coveralls badges for services no longer in use.

### Fixed
- `plot_peak_hold()` raised a `ValueError` on every call. It passed the
  `b` argument to `Axes.grid()`, which matplotlib removed in 3.5.
- `freq_bin()` returned the lower bin for a frequency falling exactly
  halfway between two bins, because the builtin `round()` rounds a
  halfway value to the nearest even bin. Halfway values now round up.
- `freq_bin()` printed debugging information on every call.
- `single_frequency_over_time()` did not subtract the frequency of the
  first bin, so it raised an `IndexError` or returned the wrong
  amplitude unless the frequency array started at 0 Hz. It now also
  rounds to the nearest bin rather than truncating, matching
  `freq_bin()`.

## v0.5.1 - 2024-12-19

### Changed
- Format the code using ruff.
- Run check and format in `just fix`.

## v0.5.0 - 2024-12-18

### Changed
- Switched to uv and Just for dependency management and task running.
- Added type hints and dropped the Python 2 `__future__` imports,
  limiting the package to Python 3.
- Build and publish using build and twine.
- Updated dependencies.

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
[ty]: https://github.com/astral-sh/ty
