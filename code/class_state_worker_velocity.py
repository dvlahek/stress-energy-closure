#!/usr/bin/env python3
"""One-shot CLASS worker with direct velocity-transfer output enabled.

This worker is deliberately short-lived.  It requests both density and velocity
transfer functions from CLASS, evaluates one PSD state, serializes only the
arrays needed by the direct-theta pair check, and exits so native CLASS memory
is returned to the OS.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from classy import Class

import hidden_channel_forensics as hf
import hidden_channel_operator_diagnostics as hd


def parse_grid(text: str):
    return [float(x) for x in text.split(",") if x.strip()]


def interp_transfer(t, key, kh):
    kin = np.asarray(t["k (h/Mpc)"], float)
    return np.interp(kh, kin, np.asarray(t[key], float))


def choose_velocity_key(keys, species: str):
    """Find the CLASS velocity-transfer column without assuming one wrapper version."""
    keys = list(keys)
    if species == "cdm":
        preferred = ("t_cdm", "theta_cdm", "v_cdm")
        tokens = ("cdm",)
    elif species == "ncdm0":
        preferred = ("t_ncdm[0]", "theta_ncdm[0]", "v_ncdm[0]")
        tokens = ("ncdm[0]", "ncdm0")
    else:
        raise ValueError(species)

    for key in preferred:
        if key in keys:
            return key

    for key in keys:
        low = key.lower().replace(" ", "")
        if any(tok in low for tok in tokens) and (
            low.startswith("t_") or low.startswith("theta_") or low.startswith("v_")
        ):
            return key
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psd", type=Path, required=True)
    ap.add_argument("--mass", type=float, required=True)
    ap.add_argument("--z-grid", required=True)
    ap.add_argument("--kh-npy", type=Path, required=True)
    ap.add_argument("--out-npz", type=Path, required=True)
    args = ap.parse_args()

    zgrid = parse_grid(args.z_grid)
    kh = np.load(args.kh_npy)

    params = hf.class_params(args.psd, args.mass, max(zgrid), max(kh), high_precision=False)
    # CLASS documents dTk and vTk as the density- and velocity-transfer outputs.
    # Keep mPk because P_cb is used in the density-weighted control.
    params["output"] = "mPk,dTk,vTk"

    c = Class()
    c.set(params)
    c.compute()
    try:
        payload = {"z_grid": np.asarray(zgrid, dtype=float)}
        transfer_keys = None
        cdm_key = None
        ncdm_key = None

        for i, z in enumerate(zgrid):
            z = float(z)
            t = c.get_transfer(z=z, output_format="class")
            if transfer_keys is None:
                transfer_keys = sorted(t.keys())
                cdm_key = choose_velocity_key(transfer_keys, "cdm")
                ncdm_key = choose_velocity_key(transfer_keys, "ncdm0")

            pk = hd.pk_cb_h(c, z, kh)
            vrel_density = hd.relative_velocity_kms(c, z, kh, 0.004)

            payload[f"pk_{i}"] = np.asarray(pk, dtype=float)
            payload[f"vrel_density_{i}"] = np.asarray(vrel_density, dtype=float)

            if cdm_key is None or ncdm_key is None:
                payload[f"theta_rel_{i}"] = np.asarray([], dtype=float)
            else:
                theta_rel = interp_transfer(t, ncdm_key, kh) - interp_transfer(t, cdm_key, kh)
                payload[f"theta_rel_{i}"] = np.asarray(theta_rel, dtype=float)

        payload["transfer_keys_json"] = np.array(json.dumps(transfer_keys or []))
        payload["cdm_velocity_key"] = np.array(cdm_key or "")
        payload["ncdm_velocity_key"] = np.array(ncdm_key or "")
        payload["theta_available"] = np.array(int(cdm_key is not None and ncdm_key is not None), dtype=np.int8)
        np.savez_compressed(args.out_npz, **payload)
    finally:
        c.struct_cleanup()
        c.empty()


if __name__ == "__main__":
    main()
