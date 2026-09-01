# Copyright (c) 2013-2026 The siganalysis developers. All rights reserved.
# Project site: https://github.com/questrail/siganalysis
# Use of this source code is governed by a MIT-style license that
# can be found in the LICENSE.txt file for the project.
"""Check a built wheel from outside the source tree.

Every other check in this project runs against `src/`, so a packaging mistake
that leaves a module or `py.typed` out of the distribution passes ruff,
pyright, and the whole suite and ships anyway. Run this with the wheel
installed somewhere `src/` cannot be reached:

    uv run --no-project --with dist/*.whl python scripts/smoke_test_wheel.py 0.8.0

Python puts this file's own directory on `sys.path` rather than the working
directory, so `import siganalysis` below can only resolve to the installed
wheel. Installing that wheel without the `plotting` extra is also the only way
to check that matplotlib is still optional: the suite runs with matplotlib
installed, so nothing there can tell whether it has become a hard dependency.
"""

import argparse
import importlib.metadata
import importlib.resources
import importlib.util

import siganalysis


def main(expected: str) -> None:
    installed = importlib.metadata.version("siganalysis")
    if installed != expected:
        raise SystemExit(f"The wheel installed {installed}, expected {expected}")

    if importlib.util.find_spec("matplotlib") is not None:
        raise SystemExit(
            "matplotlib is installed here, so this cannot tell whether it is "
            "still optional; install the wheel without the plotting extra"
        )

    # The package's own list, so that a plotting function added later is
    # checked as one rather than being expected to import without matplotlib.
    plotting_names = siganalysis._PLOTTING_NAMES

    for name in sorted(set(siganalysis.__all__) - plotting_names):
        getattr(siganalysis, name)

    # Reaching a plotting name without matplotlib has to fail, and the failure
    # has to name the extra that fixes it.
    for name in sorted(plotting_names):
        try:
            getattr(siganalysis, name)
        except ImportError as error:
            if "plotting" not in str(error):
                raise SystemExit(
                    f"{name} failed without naming the plotting extra: {error}"
                ) from error
        else:
            raise SystemExit(f"{name} resolved without matplotlib installed")

    if not importlib.resources.files("siganalysis").joinpath("py.typed").is_file():
        raise SystemExit("The wheel is missing py.typed")

    print(f"siganalysis {installed} imported from {siganalysis.__file__}")
    print("matplotlib is still optional")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check a built siganalysis wheel from outside the source tree."
    )
    parser.add_argument(
        "expected_version", help="the version the wheel is expected to install"
    )
    main(parser.parse_args().expected_version)
