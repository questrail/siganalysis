# Copyright (c) 2013-2026 The siganalysis developers. All rights reserved.
# Project site: https://github.com/questrail/siganalysis
# Use of this source code is governed by a MIT-style license that
# can be found in the LICENSE.txt file for the project.
from typing import Any

from .siganalysis import (
    STFT_WINDOWS,
    WINDOW_FUNCTIONS,
    __version__,
    calculate_peak_hold,
    freq_bin,
    hz2khz,
    single_frequency_over_time,
    smooth,
    smooth2,
    stft,
    time_slice_zip,
)

# The plotting functions are not imported here, since importing them imports
# matplotlib, which is only an optional dependency. They are reached through
# the lazy lookup below, so that they stay available as siganalysis.plot_*
# without a plain import siganalysis needing matplotlib at all.
_PLOTTING_NAMES = frozenset({"plot_peak_hold", "plot_spectrogram"})

__all__ = [
    "STFT_WINDOWS",
    "WINDOW_FUNCTIONS",
    "__version__",
    "calculate_peak_hold",
    "freq_bin",
    "hz2khz",
    "plot_peak_hold",
    "plot_spectrogram",
    "single_frequency_over_time",
    "smooth",
    "smooth2",
    "stft",
    "time_slice_zip",
]


def __getattr__(name: str) -> Any:
    """Import the plotting module on the first use of a name from it.

    Python calls this for a name the module does not already hold, so the
    cost of importing matplotlib is paid by whoever plots and by nobody else.
    """
    if name in _PLOTTING_NAMES:
        from . import plotting

        return getattr(plotting, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """List the plotting names too, which __getattr__ would otherwise hide."""
    return sorted(__all__)
