#!/usr/bin/env python3
"""Compatibility wrapper for the R cross-score visualiser.

The maintained implementation is compose_crossscore.R. Keeping this wrapper
prevents stale Python plotting code from drifting away from the paper figure
style or metric definition.
"""

from __future__ import annotations

import os
import subprocess
import sys


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    script = os.path.join(here, "compose_crossscore.R")
    cmd = ["Rscript", script, *sys.argv[1:]]
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
