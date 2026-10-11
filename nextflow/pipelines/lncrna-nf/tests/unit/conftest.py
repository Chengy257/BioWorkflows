"""Make the pipeline-level bin/ helper scripts importable for the unit tests.

The helper scripts live in <pipeline>/bin/ (module-level bin/ is not
auto-staged in this layout - see the validation notes); they are plain
scripts without a package, so the tests add the directory to sys.path.
"""
import os
import sys

BIN_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "bin"))
if BIN_DIR not in sys.path:
    sys.path.insert(0, BIN_DIR)
