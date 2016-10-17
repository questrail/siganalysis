# siganalysis

[![PyPi Version][pypi ver image]][pypi ver link]
[![Build Status][travis image]][travis link]
[![Coverage Status][coveralls image]][coveralls link]
[![License Badge][license image]][LICENSE.txt]

Python (2.7+/3.4+) routines for analyzing signals. Some of the functions
include:

- Calculating [Short-Time Fourier Transform][stft]
- Smoothing a signal
- Plotting an STFT's spectrogram
- Calculating the peak hold of an STFT in the freq domain
- Plotting the peak hold of an STFT

The above functions are handy when analyzing signals recorded in the
time domain, such as using a TEAC LX-10 data recorder, and seeing the
frequency spectrum. This is useful for Electromagnetic Compatibility
(EMC) analyses.

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

[siganalysis][] is developed using [Scott Chacon][]'s [GitHub Flow][].
To contribute, fork [siganalysis][], create a feature branch, and then
submit a pull request.  [GitHub Flow][] is summarized as:

- Anything in the `master` branch is deployable
- To work on something new, create a descriptively named branch off of
  `master` (e.g., `new-oauth2-scopes`)
- Commit to that branch locally and regularly push your work to the same
  named branch on the server
- When you need feedback or help, or you think the branch is ready for
  merging, open a [pull request][].
- After someone else has reviewed and signed off on the feature, you can
  merge it into master.
- Once it is merged and pushed to `master`, you can and *should* deploy
  immediately.

### Running the Tests

[siganalysis][] has been tested against Python 2.7 and Python 3.4. While
[siganalysis][] may work with other versions of Python, it hasn't been
tested with anything else at this time. Assuming you have [virtualenv][]
and [virtualenvwrapper][] installed, you can run the following commands
to clone the GitHub repository, create a virtualenv, and run the unit
tests:

```bash
$ git clone https://github.com/questrail/siganalysis.git
$ cd siganalysis
$ mkvirtualenv siganalysis
$ pip install -r requirements.txt
$ nosetests
```

## License

[siganalysis][] is released under the MIT license. Please see the
[LICENSE.txt][] file for more information.

[coveralls image]: http://img.shields.io/coveralls/questrail/siganalysis/master.svg
[coveralls link]: https://coveralls.io/r/questrail/siganalysis
[github flow]: http://scottchacon.com/2011/08/31/github-flow.html
[LICENSE.txt]: https://github.com/questrail/siganalysis/blob/develop/LICENSE.txt
[license image]: http://img.shields.io/pypi/l/siganalysis.svg
[numpy]: http://www.numpy.org
[matplotlib]: http://matplotlib.org
[pull request]: https://help.github.com/articles/using-pull-requests
[pypi ver image]: http://img.shields.io/pypi/v/siganalysis.svg
[pypi ver link]: https://pypi.python.org/pypi/siganalysis/
[scipy]: http://www.scipy.org
[scott chacon]: http://scottchacon.com/about.html
[siganalysis]: https://github.com/questrail/siganalysis
[stft]: http://en.wikipedia.org/wiki/Short-time_Fourier_transform
[travis image]: http://img.shields.io/travis/questrail/siganalysis/master.svg
[travis link]: https://travis-ci.org/questrail/siganalysis
[virtualenv]: https://virtualenv.pypa.io/en/latest/
[virtualenvwrapper]: http://virtualenvwrapper.readthedocs.org/en/latest/
