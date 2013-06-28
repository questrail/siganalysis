# siganalysis

Routines for analyzing signals in Python. Some of the functions include:

* Calculating Short-Term Fourier Transform
* Smoothing a signal
* Plotting an STFT's spectrogram
* Calculating the peak hold of an STFT in the freq domain
* Plotting the peak hold of an STFT

The above functions are handy when analyzing signals recorded in the
time domain, such as using a TEAC LX-10 data recorder, and seeing the
frequency spectrum. This is usefull for Electromagnetic Compatibiliity
(EMC) analyses.

## Installation

You can install [siganalysis] either via the Python Package Index (PyPI)
or from source.

To install using pip:

```bash
$ pip install siganalysis
```

**Source:** https://github.com/questrail/siganalysis

## Requirements

[siganalysis] requires the following Python packages:

* [numpy]
* [scipy]
* [matplotlib]

## Contributing

[siganalysis] is developed using [git-flow], which are "git extensions
to provide high-level repository operations for [Vincent Driessen's
branching model][nvie-git]." To contribute, [install git-flow], fork
[siganalysis], and then run:

```bash
$ git clone git@github.com/<username>/siganalysis.git
$ cd siganalysis
$ git branch master origin/master
$ git flow init -d
$ git flow feature start <your_feature>
```

When you're done coding and committing the changes for `your_feature`,
issue:

```bash
$ git flow feature publish <your_feature>
```

Then open a pull request to `your_feature` branch.

## License

[siganalysis] is released under the MIT license. Please see the
[LICENSE.txt] file for more information.

[siganalysis]: https://github.com/questrail/siganalysis
[numpy]: http://www.numpy.org
[scipy]: http://www.scipy.org
[matplotlib]: http://matplotlib.org
[LICENSE.txt]: https://github.com/questrail/siganalysis/blob/develop/LICENSE.txt
[git-flow]: https://github.com/nvie/gitflow
[nvie-git]: http://nvie.com/posts/a-successful-git-branching-model/
[install git-flow]: https://github.com/nvie/gitflow/wiki/Installation
