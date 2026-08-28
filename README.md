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

**Source:** https://github.com/questrail/siganalysis


## Requirements

[siganalysis][] requires the following Python packages:

- [numpy][]
- [scipy][]
- [matplotlib][]


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


## License

[siganalysis][] is released under the MIT license. Please see the
[LICENSE.txt][] file for more information.


[just]: https://github.com/casey/just
[LICENSE.txt]: https://github.com/questrail/siganalysis/blob/develop/LICENSE.txt
[license image]: http://img.shields.io/pypi/l/siganalysis.svg
[numpy]: http://www.numpy.org
[matplotlib]: http://matplotlib.org
[pull request]: https://help.github.com/articles/using-pull-requests
[pypi ver image]: http://img.shields.io/pypi/v/siganalysis.svg
[pypi ver link]: https://pypi.python.org/pypi/siganalysis/
[scipy]: http://www.scipy.org
[siganalysis]: https://github.com/questrail/siganalysis
[stft]: http://en.wikipedia.org/wiki/Short-time_Fourier_transform
[uv]: https://docs.astral.sh/uv/
