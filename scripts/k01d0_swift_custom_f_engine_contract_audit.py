#!/usr/bin/env python3
"""K01D0 pinned SWIFT engine-contract audit.

Reads a separately checked-out SWIFT tree at the preregistered commit and
verifies the exact source-level contract required for custom nonthermal
PartType6 neutrino ICs. No simulations or observations are used.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
PREREG = ROOT / "source_data/k01d0_swift_custom_f_engine_contract_prereg_2026-11-06.json"
SWIFT = ROOT / "external" / "SWIFT"
OUT = ROOT / "k01d0_swift_engine_contract_result.json"


def git_blob(path: Path) -> str:
    raw = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError("K01D0_STOP: " + msg)


p = json.loads(PREREG.read_text())
need(p["stage"] == "K01D0_SWIFT_CUSTOM_F_NEUTRINO_ENGINE_CONTRACT_PREREG", "stage")
need(SWIFT.is_dir(), "external SWIFT checkout missing")

head = subprocess.check_output(
    ["git", "-C", str(SWIFT), "rev-parse", "HEAD"], text=True
).strip()
need(head == p["swift"]["commit"], f"SWIFT commit mismatch {head}")

blob_checks = {}
for rel, expected in p["swift"]["pinned_git_blobs"].items():
    path = SWIFT / rel
    got = git_blob(path)
    blob_checks[rel] = {"expected": expected, "got": got, "pass": got == expected}
    need(got == expected, f"SWIFT source blob changed: {rel}")

neutrino_h = (SWIFT / "src/neutrino/Default/neutrino.h").read_text()
neutrino_c = (SWIFT / "src/neutrino/Default/neutrino.c").read_text()
neutrino_io = (SWIFT / "src/neutrino/Default/neutrino_io.h").read_text()
nu_doc = (SWIFT / "doc/RTD/source/Neutrinos/index.rst").read_text()
ic_doc = (SWIFT / "doc/RTD/source/InitialConditions/index.rst").read_text()

checks = {}

checks["PartType6_supported"] = (
    "PartType6" in nu_doc
    and "PartType6" in ic_doc
    and "swift_type_neutrino" in ic_doc
)

checks["IC_coordinates_velocities_masses_particleIDs_read"] = all(
    token in ic_doc
    for token in ("Coordinates", "Velocities", "Masses", "ParticleIDs")
) and "peculiar" in ic_doc

checks["generate_ics_zero_leaves_input_neutrino_velocity_unchanged"] = (
    "if (!e->neutrino_properties->generate_ics) return;" in neutrino_h
    and "gp->v_full[0]" in neutrino_h
)

checks["delta_f_zero_does_not_reconstruct_initial_FD_momentum"] = (
    "if (e->neutrino_properties->use_delta_f ||" in neutrino_io
    and "We don't know what the initial momentum was and we don't need it" in neutrino_io
)

# Plain particle mode is allowed because the particle-incompatible branch only
# forbids use_linear_response/use_model_none with particles, while the
# missing-particle error only applies to delta-f/generate_ics modes.
checks["plain_particle_neutrino_mode_allowed_with_particles_present"] = (
    "if ((use_df || use_df_mesh || genics) && !with_neutrinos)" in neutrino_c
    and "else if ((use_linres || use_none) && with_neutrinos)" in neutrino_c
    and "if (!use_df && !use_df_mesh) return;" in neutrino_c
)

checks["relativistic_neutrino_drift_available"] = (
    "Relativistic Drift" in nu_doc
    and "relativistic drift correction" in nu_doc
    and "a^2" in nu_doc
)

checks["custom_Fplus_Fminus_compatible_in_principle"] = all(
    checks[k]
    for k in (
        "PartType6_supported",
        "IC_coordinates_velocities_masses_particleIDs_read",
        "generate_ics_zero_leaves_input_neutrino_velocity_unchanged",
        "delta_f_zero_does_not_reconstruct_initial_FD_momentum",
        "plain_particle_neutrino_mode_allowed_with_particles_present",
        "relativistic_neutrino_drift_available",
    )
)

for key, expected in p["required_contract"].items():
    need(checks.get(key) is expected, f"contract requirement failed: {key}")

result = {
    "stage": p["stage"],
    "status": "PASS_SWIFT_CUSTOM_F_PARTICLE_NEUTRINO_ENGINE_CONTRACT",
    "swift_commit": head,
    "blob_checks": blob_checks,
    "contract_checks": checks,
    "forbidden_K01D1_modes": p["forbidden_for_K01D1"],
    "interpretation": p["interpretation"],
    "guards": p["guards"],
    "next_stage": "K01D1_EXECUTABLE_SWIFT_CUSTOM_F_COSMOLOGICAL_TRANSPORT_SMOKE",
}

OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
print("K01D0_SWIFT_COMMIT_PASS", head, flush=True)
for k, v in checks.items():
    print("K01D0_CONTRACT", k, v, flush=True)
print("K01D0_SWIFT_CUSTOM_F_ENGINE_CONTRACT_PASS", flush=True)
print("K01D0_RESULT", OUT, flush=True)
