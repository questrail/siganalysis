from distutils.core import setup
setup(
    name='siganalysis',
    version='0.1.0',
    author='Matthew Rankin',
    author_email='matthew@questrail.com',
    py_modules=['siganalysis'],
    url='http://github.com/questrail/siganalysis',
    license='LICENSE.txt',
    description='Perform signal analysis',
    requires=['numpy (>=1.6.0)',
              'scipy (>=0.11.0)'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Development Status :: 3 - Alpha',
        'Operating System :: OS Independent',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
