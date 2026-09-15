#!/usr/bin/env python3
"""Low-memory launcher for hidden_channel_forensics.py.

The scientific calculation is unchanged.  The only difference is execution
isolation: every CLASS state is computed in a fresh Python subprocess via
class_state_worker.py.  This prevents native memory retained by repeated classy
calls from accumulating in the long-lived parent process, which can otherwise
trigger the Linux/WSL OOM killer after many CLASS evaluations.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

import hidden_channel_forensics as hf

WORKER = Path(__file__).with_name("class_state_worker.py")


def isolated_build_state(psd: Path, mass: float, zgrid, kh, high_precision=False):
    psd = Path(psd)
    workdir = psd.parent
    workdir.mkdir(parents=True, exist_ok=True)

    kh_fd, kh_name = tempfile.mkstemp(prefix="hf_kh_", suffix=".npy", dir=workdir)
    os.close(kh_fd)
    out_fd, out_name = tempfile.mkstemp(prefix="hf_state_", suffix=".npz", dir=workdir)
    os.close(out_fd)
    kh_path = Path(kh_name)
    out_path = Path(out_name)

    try:
        np.save(kh_path, np.asarray(kh, dtype=float))
        cmd = [
            sys.executable,
            str(WORKER),
            "--psd", str(psd.resolve()),
            "--mass", repr(float(mass)),
            "--z-grid", ",".join(repr(float(z)) for z in zgrid),
            "--kh-npy", str(kh_path.resolve()),
            "--out-npz", str(out_path.resolve()),
        ]
        if high_precision:
            cmd.append("--high-precision")
        subprocess.run(cmd, check=True)

        with np.load(out_path, allow_pickle=False) as d:
            theta_available = bool(int(d["theta_available"]))
            transfer_keys = json.loads(str(d["transfer_keys_json"].item()))
            zs = np.asarray(d["z_grid"], dtype=float)
            state = {}
            for i, z in enumerate(zs):
                theta = np.asarray(d[f"theta_{i}"], dtype=float)
                state[float(z)] = {
                    "vrel": np.asarray(d[f"vrel_{i}"], dtype=float).copy(),
                    "pk": np.asarray(d[f"pk_{i}"], dtype=float).copy(),
                    "theta_rel": None if theta.size == 0 else theta.copy(),
                }
        return state, theta_available, transfer_keys
    finally:
        for p in (kh_path, out_path):
            try:
                p.unlink()
            except FileNotFoundError:
                pass


def main():
    # Monkey-patch only the state-construction routine.  All operators,
    # normalization, validation and summary logic remain exactly the same as in
    # hidden_channel_forensics.py.
    hf.build_state = isolated_build_state
    hf.main()


if __name__ == "__main__":
    main()
