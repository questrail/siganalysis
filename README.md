# siganalysis

[![PyPI Version][pypi ver image]][pypi ver link]
[![Python Versions][pyversions image]][pypi ver link]
[![CI Status][ci image]][ci link]
[![Coverage Status][coveralls image]][coveralls link]
[![License Badge][license image]][LICENSE.txt]

Python (3.12+) routines for analyzing signals. Some of the functions include:

- Calculating [Short-Time Fourier Transform][stft]
- Smoothing a signal
- Plotting an STFT's spectrogram
- Calculating the peak hold of an STFT in the freq domain
- Plotting the peak hold of an STFT

The above functions are handy when analyzing signals recorded in the time
domain, such as using a TEAC LX-10 data recorder, and seeing the frequency
spectrum for Electromagnetic Compatibility (EMC) analyses.


## Installation

You can install [siganalysis][] either via the Python Package Index
(PyPI) or from source.

To add it to a project managed with [uv][], which records it in your
`pyproject.toml` and lock file:

```bash
$ uv add siganalysis
```

Or to install it with pip:

```bash
$ pip install siganalysis
```

The plotting functions need [matplotlib][], which is an optional
dependency so that it is only installed for those who plot. Install it
alongside [siganalysis][] with the `plotting` extra:

```bash
$ uv add "siganalysis[plotting]"
$ pip install siganalysis[plotting]
```

**Source:** https://github.com/questrail/siganalysis


## Requirements

[siganalysis][] requires the following Python packages:

- [numpy][]
- [scipy][]

The plotting functions, `plot_spectrogram()` and `plot_peak_hold()`,
additionally require [matplotlib][], which is installed with the
`plotting` extra described above. They live in `siganalysis.plotting`,
so that importing [siganalysis][] does not import [matplotlib][], but
they remain reachable from the package itself:

```python
import siganalysis

siganalysis.plot_peak_hold(axis, stft_data, freq_array)
```


## Release Notes

The [CHANGELOG][] records every change, and each release posts its own
section of it to the [releases page][releases]. That posting is done by
the release workflow rather than by hand, so the CHANGELOG is the one
place a change has to be written down. The longer notes kept in
[docs/releases][release notes] cover v0.8.0, which predates that.


## Contributing

Contributions are welcome! To contribute please:

1. Fork the repository
2. Create a feature branch
3. Add code and tests
4. Pass lint and tests
5. Submit a [pull request][]


## Development Setup

[siganalysis][] uses [uv][] to manage the virtualenv and dependencies,
and [just][] as the task runner.

```bash
$ brew install uv just
```

With [uv][] and [just][] installed, `uv sync` creates the virtualenv and
installs the dependencies, including the development group, and `just`
on its own lists the available recipes.

```bash
$ uv sync
$ just
```

The most common recipes are:

```bash
$ just test    # Run the tests using pytest
$ just lint    # Check lint, formatting, types, and workflows
$ just fix     # Lint and format the code using ruff, applying fixes
$ just cov     # Run the tests and report coverage
$ just add X   # Add X as a dependency
$ just out     # List the outdated dependencies
```

[ruff][] and [pyright][] are deliberately absent from that `brew
install` line. Both are dev dependencies pinned in `uv.lock` and reached
through `uv run`, so every recipe and every CI job uses the same
version. A `brew install ruff` would put a second, unpinned copy on the
path for an editor to find, and ruff releases change how code is
formatted: the editor would then reformat code that `ruff format
--check` rejects on the next run.


## Making a Release

`just release` cuts the release. It first checks that a release is
possible at all, then lints, type checks, and tests, then shows the
entries waiting under Unreleased and the version each kind of bump would
produce, and asks which to cut. Once answered it bumps the version,
closes out the CHANGELOG, updates the lock file, commits, and tags.
Pushing the tag is what publishes.

```bash
$ just release

Releasing from 0.8.0, with these entries under Unreleased:

    ### Added

    - `stft()` accepts a choice of window through its new `window` argument.

    1) patch   0.8.0 -> 0.8.1
    2) minor   0.8.0 -> 0.9.0
    3) major   0.8.0 -> 1.0.0
    q) cancel

Which release? [1] 2

Tagged v0.9.0. Publish it with:

    git push --follow-tags
```

The entries decide the bump, so the prompt puts them next to the
versions they would produce rather than leaving the choice to memory.
Answering `q`, or anything unrecognized, changes nothing.

The tag push runs the [release workflow][], which waits on the whole [CI
workflow][ci link] before it does anything else: the 3.12, 3.13, and
3.14 matrix and the dependency floor job. `git push --follow-tags` starts
both at once, so without that wait an upload could go out while a leg of
the matrix was still running, or already red. It then checks that the
tagged commit is on `master`, since a tag is only a pointer and one
placed anywhere else would otherwise publish whatever it points at,
rechecks the tag against the version in `pyproject.toml`, and builds.

Every check to that point runs against the source tree with matplotlib
installed, so the workflow then installs the wheel it just built, without
the `plotting` extra, somewhere `src/` is not on the path. That is the
only step that can catch a packaging mistake which left something out of
the distribution, and the only one that can tell whether matplotlib is
still optional. It uploads once that passes. There is no PyPI API token
anywhere: the workflow authenticates with [trusted publishing][], which
mints a short lived credential from the GitHub OIDC identity of that
run. That same identity signs a [PEP 740][] attestation for each
distribution, which PyPI serves beside the file it attests: trusted
publishing establishes who uploaded, and the attestation establishes
what was uploaded and which workflow built it. The upload skips anything
PyPI already holds, so a run that uploaded one distribution and then
failed on the other can be retried instead of stranding a version number
that PyPI will never allow to be reused.

Uploading is followed by a [GitHub release][releases] for the tag,
carrying the CHANGELOG section for that version as its notes and the
built distributions as its assets. The notes are collected before the
upload rather than after, so that a CHANGELOG with no section for the
version being released stops the release while stopping it is still
possible.

Pushing the tag is the point of no return, since PyPI never lets a
version number be reused. Everything `just release` does is local and
amendable until then, and it refuses to start against a dirty working
tree, off `master`, on a `master` behind its upstream, with a CHANGELOG
whose Unreleased section is empty, or when the tag it would create
already exists. Those refusals come before the lint and test run, so a
release that cannot happen is turned away at once rather than after the
suite. `just release-check` runs them on their own, and a refusal leaves
the version and the CHANGELOG untouched.

`just build` runs the same checks and produces the same distributions
without releasing anything, which is the way to inspect what CI would
upload.

This depends on one piece of configuration that lives outside the
repository. A [trusted publisher][trusted publishing] has to be
registered for `siganalysis` on PyPI, pointing at the
`questrail/siganalysis` repository, the `release.yml` workflow, and the
`pypi` environment. It is a one time setup per project.


## License

[siganalysis][] is released under the MIT license. Please see the
[LICENSE.txt][] file for more information.


[CHANGELOG]: https://github.com/questrail/siganalysis/blob/master/CHANGELOG.md
[ci image]: https://github.com/questrail/siganalysis/actions/workflows/ci.yml/badge.svg?branch=master
[ci link]: https://github.com/questrail/siganalysis/actions/workflows/ci.yml
[coveralls image]: https://coveralls.io/repos/github/questrail/siganalysis/badge.svg?branch=master
[coveralls link]: https://coveralls.io/github/questrail/siganalysis?branch=master
[just]: https://just.systems/
[license image]: https://img.shields.io/pypi/l/siganalysis.svg
[LICENSE.txt]: https://github.com/questrail/siganalysis/blob/master/LICENSE.txt
[matplotlib]: http://matplotlib.org
[numpy]: http://www.numpy.org
[PEP 740]: https://peps.python.org/pep-0740/
[pull request]: https://help.github.com/articles/using-pull-requests
[pypi ver image]: https://img.shields.io/pypi/v/siganalysis.svg
[pypi ver link]: https://pypi.python.org/pypi/siganalysis/
[pyright]: https://microsoft.github.io/pyright/
[pyversions image]: https://img.shields.io/pypi/pyversions/siganalysis.svg
[release notes]: https://github.com/questrail/siganalysis/tree/master/docs/releases
[release workflow]: https://github.com/questrail/siganalysis/blob/master/.github/workflows/release.yml
[releases]: https://github.com/questrail/siganalysis/releases
[ruff]: https://docs.astral.sh/ruff/
[scipy]: http://www.scipy.org
[siganalysis]: https://github.com/questrail/siganalysis
[stft]: http://en.wikipedia.org/wiki/Short-time_Fourier_transform
[trusted publishing]: https://docs.pypi.org/trusted-publishers/
[uv]: https://docs.astral.sh/uv/
