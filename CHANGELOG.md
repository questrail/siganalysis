# CHANGELOG.md
This file contains all notable changes to the [siganalysis][] project.

## Unreleased

### Added

- Continuous integration on GitHub Actions, which this project had none
  of. Every push and pull request now lints, checks formatting, type
  checks, and runs the suite on 3.12 and 3.13, the versions the
  `pyproject.toml` classifiers claim. A second job installs the oldest
  numpy, scipy, and matplotlib that `pyproject.toml` allows, with
  `--resolution lowest-direct`, so that the version floors are a
  tested promise rather than a hopeful one; the lock file pins the
  newest of each, so nothing else exercises them. A third audits the
  workflows with [zizmor][], the part of the repository that can mint a
  PyPI credential having otherwise been read by eye alone. That job
  restores the uv cache without saving it, since it installs what the
  3.13 leg installs and the two would otherwise race to write one key.
  Coverage goes to Coveralls from the 3.13 leg.
- Releases publish from a tag rather than from a laptop. `just release`
  refuses a dirty tree, a branch other than `master`, a `master` behind
  its upstream, an empty Unreleased section, or an existing tag; then
  lints and tests; then shows the entries waiting to ship beside the
  version each kind of bump would produce, and asks which to cut. It
  bumps the version, closes out the CHANGELOG, commits, and tags.
  Pushing the tag publishes. `just release-check` runs the refusals on
  their own.
- The release workflow waits on the whole CI run before it uploads
  anything, confirms the tag sits on `master` and matches the version
  in `pyproject.toml`, and authenticates to PyPI with
  [trusted publishing][], so there is no API token to paste, store, or
  leak. It signs a [PEP 740][] attestation for each distribution
  against the same identity, and creates a GitHub release carrying the
  CHANGELOG section for that version as its notes.
- Dependabot keeps the pinned actions and the lock file moving. The
  actions in both workflows are pinned to commit SHAs, which a fix
  published upstream does not reach on its own the way a moving tag
  would, so pinning without something to move it would amount to
  staying on one commit forever. Dependabot reads `pyproject.toml` and
  `uv.lock` together as well, so a dependency update arrives as a lock
  file change that CI checks with `uv sync --locked` rather than as a
  resolution done on the runner.
- `scripts/smoke_test_wheel.py`, which installs the built wheel where
  `src/` cannot be reached and without the `plotting` extra, then
  checks the version, every public name, and `py.typed`. Every other
  check runs against the source tree with matplotlib installed, so this
  is the only one that can catch a packaging mistake, and the only one
  that can tell whether matplotlib is still optional. `just build` and
  the release workflow run the same command.
- `py.typed`, so that the type hints already written reach anyone
  installing the package rather than stopping at this repository.
- A coverage floor of 95%, which is what the suite covers today. Below
  that the run fails, so uncovered code has to arrive with either a
  test or a deliberate edit to the floor.
- `just doc`, `just up-all`, and `just lint`, and a `Justfile` that is
  otherwise line for line the one in [applyaf][], so that moving
  between the two projects does not mean learning a second set of
  recipe names.

### Changed

- Type checking is done by [pyright][] rather than [ty][], which is
  still a 0.0.x release, and it runs inside `just lint` rather than as
  a separate `just check` that every caller had to remember. Reaching
  it through `lint` means `just build` and `just release` cannot skip
  it. [scipy-stubs][] is a new dev dependency: without it pyright
  cannot see the window functions in `scipy.signal.windows`.
- `plot_spectrogram()` and `plot_peak_hold()` are declared to a type
  checker under `TYPE_CHECKING`, so an editor offers their real
  signatures instead of names that appear from nowhere. They still
  arrive through the module's `__getattr__` at run time, so a plain
  `import siganalysis` still does not import matplotlib.
- The ruff rule set is selected explicitly in `pyproject.toml` rather
  than left at the default, so that it is a deliberate choice. The
  handful of places the new rules object to are deliberate, and each
  carries a `noqa` saying so.
- The license is declared as an SPDX expression with `license-files`,
  which is what replaced the `License ::` classifier that used to
  carry it.
- `.python-version` is tracked rather than ignored. It decides the
  interpreter a contributor's `uv sync` builds against, and the file
  was being read while being excluded from the repository.

### Removed

- `just deploy`, which published from a laptop against a pasted PyPI
  token, and `just test-all`, which ran the suite against each
  supported Python locally because there was no CI to do it. The
  release workflow and the CI matrix replace them.
- `AUTHORS.md`, along with the pointer to it in the copyright notice.
  The repository history is the record of who wrote what.

## v0.8.0 - 2026-08-31

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
- The plotting functions live in `siganalysis.plotting`, and matplotlib
  is now an optional dependency installed with the `plotting` extra, so
  that it is only installed for those who plot. Closes #5.

  `siganalysis.plot_spectrogram()` and `siganalysis.plot_peak_hold()`
  still resolve, so no caller has to change. The package imports the
  plotting module on the first use of either name, which keeps
  matplotlib out of a plain `import siganalysis` while leaving both
  reachable where they always were. Plotting without matplotlib
  installed raises an ImportError naming the extra to install.
- The functions taking STFT data now check it against the vectors
  describing its axes, and say what does not line up rather than failing
  further down with an error about numpy indexing. `plot_spectrogram()`
  had accepted a time or frequency vector that did not match the data at
  all, and plotted it against mislabelled axes. `smooth2()` gained the
  checks `smooth()` already had, having quietly returned an empty array
  for a window longer than the signal. `plot_peak_hold()` checks that a
  `limit_array` carries the fields it is read by. Closes #12.
- `time_slice_zip()` takes an optional `minimum_samples_in_last_slice`.
  A last slice shorter than that is folded into the slice before it, so
  that a sample count just past a multiple of the slice size no longer
  leaves a last slice too short to process. The samples are kept rather
  than dropped, so the last slice grows instead. Closes #20.
- `time_slice_zip()` rejects fewer than one sample per time slice, which
  had put it into an endless loop.
- The examples in the docstrings now run as doctests under pytest, so an
  example that stops working fails the test suite. The example in
  `smooth()` had been broken for some time: `np.linspace(-2, 2, 0.1)`
  raises a TypeError, since the third argument is a count rather than a
  step. Closes #10.
- Tests for `smooth()`, `smooth2()`, and `plot_peak_hold()`, the three
  functions that still had none. Closes #4. The suite has grown from 53
  tests to 154.
- `WINDOW_FUNCTIONS`, the windows `smooth()` accepts, is now exported
  alongside `STFT_WINDOWS`.
- `siganalysis.__version__` reports the installed version. It was only
  reachable as `siganalysis.siganalysis.__version__` before.

### Changed
- **`stft()` takes `window` in place of `use_hamming_window`.** The
  argument is the name of a window or None for no window, so
  `use_hamming_window=False` becomes `window=None` and
  `use_hamming_window=True` becomes `window="hamming"`, which is also
  the default. Passing the old argument raises a TypeError.
- **`stft()` returns different amplitudes when a window is used.** See
  the fix below. Amplitudes are about 7% (0.594 dB) lower than they
  were, and the amplitude reported for a tone is now the same with a
  window on as with it off.
- **matplotlib is no longer installed by default.** An existing install
  that plots needs `pip install siganalysis[plotting]` on upgrade.
- **`plot_spectrogram()` draws a slightly different range** for a plot
  range that does not land on a bin. The bin holding a value is now the
  one whose value is nearest to it, rather than the one found by
  truncating. Closes #11.

### Fixed
- `stft()` scaled a Hamming windowed frame by 2 rather than by the gain
  of the window. Applying a window scales every amplitude by the mean of
  its samples, its coherent gain, which is 0.5354 for a Hamming window,
  so the correction needed is 1/0.5354 = 1.8678. Scaling by 2 instead
  left every windowed amplitude 7.08% (0.594 dB) high, which is a
  meaningful error when comparing a peak hold against a limit line. A
  1.0 amplitude tone on the center of a bin now reads 1.0 rather than
  1.0708.
- `_bin_holding()` truncated while `freq_bin()` rounded to nearest, so
  the module held two different answers to which bin holds a value. For
  10 Hz bins, 6 Hz gave the 0 Hz bin one way and the 10 Hz bin the
  other. Both round to nearest now, which is also exactly the set of
  bins overlapping a requested plot range.
- `smooth()` returned a signal one sample short of the one given when
  the window length was even and equal to the length of the signal. The
  window is made odd before the length is checked now, so such a call
  raises IndexError rather than quietly returning the wrong length.

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
[applyaf]: https://github.com/questrail/applyaf
[PEP 740]: https://peps.python.org/pep-0740/
[pyright]: https://microsoft.github.io/pyright/
[scipy-stubs]: https://github.com/scipy/scipy-stubs
[siganalysis]: https://github.com/questrail/siganalysis
[trusted publishing]: https://docs.pypi.org/trusted-publishers/
[ty]: https://github.com/astral-sh/ty
[zizmor]: https://docs.zizmor.sh/
