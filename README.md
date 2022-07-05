# siganalysis
[pypi ver image]: http://img.shields.io/pypi/v/applyaf.svg
[pypi ver link]: https://pypi.python.org/pypi/applyaf

[![PyPi Version][pypi ver image]][pypi ver link]
[![Build Status][travis image]][travis link]
[![Coverage Status][coveralls image]][coveralls link]
[![License Badge][license image]][LICENSE.txt]

Python (3.6+) routines for analyzing signals. Some of the functions include:

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

### Development Setup Using pyenv

Use the following commands to create a Python 3.9.9 virtualenv using [pyenv][]
and [pyenv-virtualenv][], install the requirements in the virtualenv named
`siganalysis`, and list the available [Invoke][] tasks.

```bash
$ pyenv virtualenv 3.9.9 siganalysis
$ pyenv activate siganalysis
$ pip install -r requirements.txt
$ inv -l
```


## License

[siganalysis][] is released under the MIT license. Please see the
[LICENSE.txt][] file for more information.


[coveralls image]: http://img.shields.io/coveralls/questrail/siganalysis/master.svg
[coveralls link]: https://coveralls.io/r/questrail/siganalysis
[invoke]: https://www.pyinvoke.org/
[LICENSE.txt]: https://github.com/questrail/siganalysis/blob/develop/LICENSE.txt
[license image]: http://img.shields.io/pypi/l/siganalysis.svg
[numpy]: http://www.numpy.org
[matplotlib]: http://matplotlib.org
[pull request]: https://help.github.com/articles/using-pull-requests
[pyenv]: https://github.com/pyenv/pyenv
[pyenv-virtualenv]: https://github.com/pyenv/pyenv-virtualenv
[pypi ver image]: http://img.shields.io/pypi/v/siganalysis.svg
[pypi ver link]: https://pypi.python.org/pypi/siganalysis/
[scipy]: http://www.scipy.org
[siganalysis]: https://github.com/questrail/siganalysis
[stft]: http://en.wikipedia.org/wiki/Short-time_Fourier_transform
[travis image]: http://img.shields.io/travis/questrail/siganalysis/master.svg
[travis link]: https://travis-ci.org/questrail/siganalysis
[virtualenv]: https://virtualenv.pypa.io/en/latest/
[virtualenvwrapper]: http://virtualenvwrapper.readthedocs.org/en/latest/
