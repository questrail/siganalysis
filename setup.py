import codecs
import os
import re

from setuptools import setup

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


try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except(IOError, ImportError):
    long_description = open('README.md').read()

setup(
    name='siganalysis',
    version=find_version('siganalysis.py'),
    author='Matthew Rankin',
    author_email='matthew@questrail.com',
    py_modules=['siganalysis'],
    url='http://github.com/questrail/siganalysis',
    license='MIT',
    description='Perform signal analysis',
    long_description=long_description,
    requires=['matplotlib (>=1.5.0)',
              'numpy (>=1.11.0)',
              'scipy (>=0.18.0)'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'License :: OSI Approved :: MIT License',
        'Development Status :: 3 - Alpha',
        'Operating System :: OS Independent',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
