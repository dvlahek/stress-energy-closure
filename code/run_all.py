#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent

commands = [
    [sys.executable, HERE / "direct_vs_memory.py", "--full"],
    [sys.executable, HERE / "mass_sweep.py"],
    [sys.executable, HERE / "hierarchy_test.py"],
]

for command in commands:
    print("\n" + "=" * 80)
    print("RUNNING:", " ".join(map(str, command)))
    subprocess.check_call(command)

print("\nAll numerical campaigns finished.")
