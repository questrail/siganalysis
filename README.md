# siganalysis

[![PyPi Version][pypi ver image]][pypi ver link]
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

To install using pip:

```bash
$ pip install siganalysis
```

The plotting functions need [matplotlib][], which is an optional
dependency so that it is only installed for those who plot. Install it
alongside [siganalysis][] with the `plotting` extra:

```bash
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

Release notes describing what changed and what a caller has to update
are posted to the [releases page][releases] and kept in
[docs/releases][release notes]. The [CHANGELOG][] records every change,
including the releases that predate those notes.


## Contributing

Contributions are welcome! To contribute please:

1. Fork the repository
2. Create a feature branch
3. Add code and tests
4. Pass lint and tests
5. Submit a [pull request][]


## Development Setup

[siganalysis][] uses [uv][] to manage the virtualenv and dependencies, and
[just][] as the task runner.

Use the following commands to create the virtualenv, install the dependencies
(including the development group), and list the available [just][] recipes.

```bash
$ uv sync
$ just
```

The most common recipes are:

```bash
$ just test    # Run the tests using pytest
$ just fix     # Lint and format the code using ruff
$ just add X   # Add X as a dependency
$ just out     # List the outdated dependencies
```


## Making a Release

Releasing is a single recipe, `just deploy`, which refuses to start
unless everything it needs is in place. Prepare the release first:

1. Bump the version with `uv version --bump minor`, which updates both
   `pyproject.toml` and `uv.lock`.
2. Add the entry for the new version to the [CHANGELOG][].
3. Write the release notes as `docs/releases/vX.Y.Z.md`.
4. Commit and push, since a release is not allowed to get ahead of the
   repository.

```bash
$ uv version --bump minor
$ git commit -am "Release vX.Y.Z"
$ git push origin master
$ just deploy
```

`just deploy` confirms the tree is ready, lints, type checks, runs the
tests against every supported Python, builds from a cleared `dist/`,
installs the built wheel on its own to confirm that it imports, and then
asks for a PyPI token and publishes. The tag is written and pushed only
after PyPI accepts the upload, so that a tag is evidence a release
shipped rather than evidence that one was attempted. Do not tag by hand.

Run `just release-check` on its own at any point to see which of those
conditions is not yet met.

Finally, post the notes to the [releases page][releases]. The first line
of the notes file is a heading that GitHub renders from the title
already, so leave it out:

```bash
$ tail -n +3 docs/releases/vX.Y.Z.md |
    gh release create vX.Y.Z --title vX.Y.Z --notes-file -
```


## License

[siganalysis][] is released under the MIT license. Please see the
[LICENSE.txt][] file for more information.


[CHANGELOG]: https://github.com/questrail/siganalysis/blob/master/CHANGELOG.md
[just]: https://github.com/casey/just
[LICENSE.txt]: https://github.com/questrail/siganalysis/blob/master/LICENSE.txt
[license image]: http://img.shields.io/pypi/l/siganalysis.svg
[numpy]: http://www.numpy.org
[matplotlib]: http://matplotlib.org
[pull request]: https://help.github.com/articles/using-pull-requests
[pypi ver image]: http://img.shields.io/pypi/v/siganalysis.svg
[pypi ver link]: https://pypi.python.org/pypi/siganalysis/
[release notes]: https://github.com/questrail/siganalysis/tree/master/docs/releases
[releases]: https://github.com/questrail/siganalysis/releases
[scipy]: http://www.scipy.org
[siganalysis]: https://github.com/questrail/siganalysis
[stft]: http://en.wikipedia.org/wiki/Short-time_Fourier_transform
[uv]: https://docs.astral.sh/uv/
