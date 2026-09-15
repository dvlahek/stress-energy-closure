#!/usr/bin/env python3
"""One-shot CLASS state worker used by hidden_channel_forensics_lowmem.py.

Each invocation performs exactly one CLASS compute, serializes the compact arrays
needed by the parent process, and exits.  Exiting the worker guarantees that any
native CLASS allocations are returned to the OS even if classy/CLASS retains
memory across repeated compute/cleanup cycles in a long-lived Python process.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import hidden_channel_forensics as hf


def parse_grid(text: str):
    return [float(x) for x in text.split(",") if x.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psd", type=Path, required=True)
    ap.add_argument("--mass", type=float, required=True)
    ap.add_argument("--z-grid", required=True)
    ap.add_argument("--kh-npy", type=Path, required=True)
    ap.add_argument("--out-npz", type=Path, required=True)
    ap.add_argument("--high-precision", action="store_true")
    args = ap.parse_args()

    zgrid = parse_grid(args.z_grid)
    kh = np.load(args.kh_npy)
    state, theta_available, transfer_keys = hf.build_state(
        args.psd, args.mass, zgrid, kh, high_precision=args.high_precision
    )

    payload = {
        "theta_available": np.array(int(theta_available), dtype=np.int8),
        "transfer_keys_json": np.array(json.dumps(transfer_keys or [])),
        "z_grid": np.asarray(zgrid, dtype=float),
    }
    for i, z in enumerate(zgrid):
        row = state[float(z)]
        payload[f"vrel_{i}"] = np.asarray(row["vrel"], dtype=float)
        payload[f"pk_{i}"] = np.asarray(row["pk"], dtype=float)
        if row["theta_rel"] is None:
            payload[f"theta_{i}"] = np.asarray([], dtype=float)
        else:
            payload[f"theta_{i}"] = np.asarray(row["theta_rel"], dtype=float)

    np.savez_compressed(args.out_npz, **payload)


if __name__ == "__main__":
    main()
