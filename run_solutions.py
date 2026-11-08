#!/usr/bin/env python3
"""Run the full ischaemic-heart-disease solution simulation from any directory."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from ihd_solutions.cli import main  # noqa: E402

if __name__ == "__main__":
    argv = sys.argv[1:] or ["simulate"]
    raise SystemExit(main(argv))
