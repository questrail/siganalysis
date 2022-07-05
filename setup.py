import codecs
import os
import re

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

here = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    """Read parts of a file

    Taken from pip's setup.py
    intentionally *not* adding an encoding option to open
    see: https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    """
    return codecs.open(os.path.join(here, *parts), 'r').read()


def find_version(*file_paths):
    """Find version in source file

    Read the version number from a source file.
    Code taken from pip's setup.py
    """
    version_file = read(*file_paths)
    # The version line must have the form:
    # __version__ = 'ver'
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setuptools.setup(
    name='siganalysis',
    version=find_version('siganalysis.py'),
    author='Matthew Rankin',
    author_email='matthew@questrail.com',
    py_modules=['siganalysis'],
    url='http://github.com/questrail/siganalysis',
    license='MIT',
    description='Perform signal analysis',
    long_description=long_description,
    long_description_content_type="text/markdown",
    requires=['matplotlib (>=3.5.1)',
              'numpy (>=1.22.3)',
              'scipy (>=1.8.0)'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Development Status :: 4 - Beta',
        'Operating System :: OS Independent',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
