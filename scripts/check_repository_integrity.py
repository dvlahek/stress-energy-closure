#!/usr/bin/env python3
"""Check the published repository without downloading survey catalogues.

This check validates source syntax, local file references, retained checksums,
and agreement among the reported observational summaries. It does not rerun
survey pair counts, mock generation, CLASS or the Einstein–Vlasov calculations.
"""
from __future__ import annotations

import ast
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
LOCAL_PREFIXES = ("code/", "scripts/", "source_data/", "docs/", "observational_forecast/", ".github/workflows/")
MARKDOWN_LINK = re.compile(r"\]\(([^)]+)\)")
BACKTICK_PATH = re.compile(r"`((?:code|scripts|source_data|docs|observational_forecast|\.github/workflows)/[\w./-]+)`")
RELATIVE_BACKTICK = re.compile(r"`((?:\.\./|\./)+[\w./-]+)`")
RUN_PATH = re.compile(r"\b(?:python(?:3)?|bash)\s+((?:code|scripts)/[\w.-]+)")


def tracked_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "-z"], cwd=ROOT, check=True, capture_output=True
    )
    return [ROOT / p.decode("utf-8") for p in result.stdout.split(b"\0") if p]


def require_file(target: Path, source: Path, problems: list[str]) -> None:
    if not target.is_file() and not target.is_dir():
        problems.append(f"{source.relative_to(ROOT)}: missing local reference {target}")


def validate_sources(files: list[Path], problems: list[str]) -> tuple[int, int, int]:
    counts = [0, 0, 0]
    for path in files:
        if path.suffix == ".py":
            counts[0] += 1
            try:
                ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            except (SyntaxError, UnicodeError) as exc:
                problems.append(f"{path.relative_to(ROOT)}: Python syntax: {exc}")
        elif path.suffix == ".json":
            counts[1] += 1
            try:
                json.loads(path.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, UnicodeError) as exc:
                problems.append(f"{path.relative_to(ROOT)}: JSON syntax: {exc}")
        elif path.suffix == ".sh":
            counts[2] += 1
            test = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
            if test.returncode:
                problems.append(f"{path.relative_to(ROOT)}: shell syntax: {test.stderr.strip()}")
    return tuple(counts)


def validate_references(files: list[Path], problems: list[str]) -> int:
    checked = 0
    for path in files:
        if path.suffix == ".md":
            content = path.read_text(encoding="utf-8")
            for match in MARKDOWN_LINK.finditer(content):
                link = unquote(match.group(1).split("#", 1)[0])
                if not link or re.match(r"^[a-z][\w+.-]*:", link, re.I):
                    continue
                require_file((path.parent / link).resolve(), path, problems)
                checked += 1
            for match in BACKTICK_PATH.finditer(content):
                require_file(ROOT / match.group(1), path, problems)
                checked += 1
            for match in RELATIVE_BACKTICK.finditer(content):
                target = match.group(1)
                if target.endswith((".md", ".json", ".py", ".csv", ".sh", ".yml", ".txt")):
                    require_file((path.parent / target).resolve(), path, problems)
                    checked += 1
        elif path.suffix == ".json":
            obj = json.loads(path.read_text(encoding="utf-8"))
            def walk(value: object) -> None:
                nonlocal checked
                if isinstance(value, str):
                    if value.startswith(LOCAL_PREFIXES) and not any(c in value for c in " <>*{}"):
                        require_file(ROOT / value, path, problems)
                        checked += 1
                elif isinstance(value, dict):
                    for child in value.values():
                        walk(child)
                elif isinstance(value, list):
                    for child in value:
                        walk(child)
            walk(obj)
        elif path.suffix == ".yml":
            for match in RUN_PATH.finditer(path.read_text(encoding="utf-8")):
                require_file(ROOT / match.group(1), path, problems)
                checked += 1
    return checked


def validate_checksums(files: list[Path], problems: list[str]) -> int:
    n = 0
    for manifest in files:
        if manifest.name != "SHA256SUMS.txt":
            continue
        for row in manifest.read_text(encoding="utf-8").splitlines():
            if not row.strip() or row.lstrip().startswith("#"):
                continue
            match = re.fullmatch(r"([a-f0-9]{64})\s+\*?(.+)", row)
            if not match:
                problems.append(f"{manifest.relative_to(ROOT)}: malformed checksum entry: {row}")
                continue
            expected, name = match.groups()
            target = (manifest.parent / name).resolve()
            if not target.is_relative_to(manifest.parent.resolve()) or not target.is_file():
                problems.append(f"{manifest.relative_to(ROOT)}: checksum target missing: {name}")
                continue
            actual = hashlib.sha256(target.read_bytes()).hexdigest()
            n += 1
            if actual != expected:
                problems.append(f"{manifest.relative_to(ROOT)}: checksum mismatch: {name}")
    return n


def same(label: str, left: float, right: float, problems: list[str]) -> None:
    if not math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-13):
        problems.append(f"reported metadata differ: {label}: {left} != {right}")


def validate_reported_results(problems: list[str]) -> None:
    data = ROOT / "source_data"
    required = (
        "phase7_validation_manifest.json",
        "phase7_desi_fullsample_perm256_summary.json",
        "lrg_elg_exact_final_manifest_2026-09-23.json",
        "lrg_elg_windowed_finite_mock_r4_120_final_2026-09-23.json",
    )
    if not all((data / p).is_file() for p in required):
        problems.append("missing primary observational summary or final manifest")
        return
    phase = json.loads((data / required[0]).read_text(encoding="utf-8"))
    phase_result = json.loads((data / required[1]).read_text(encoding="utf-8"))
    src = phase["primary_observational_inference"]
    fit = phase_result["conservative_fit"]
    for key, name in (("conservative_wake_amplitude", "wake_amplitude"),
                      ("conservative_wake_sigma", "wake_sigma"),
                      ("conservative_wake_z", "wake_z"),
                      ("empirical_two_sided_pvalue", "empirical_two_sided_pvalue")):
        same("phase7 " + key, src[key], fit[name], problems)
    same("phase7 global permutation p", src["global_permutation_pvalue"],
         phase_result["permutation_control"]["global_empirical_pvalue"], problems)
    if src["permutation_count"] != phase_result["permutation_control"]["count"]:
        problems.append("phase7 permutation count differs between manifest and result")
    if phase["blinding"].get("fully_blinded") is not False:
        problems.append("the phase7 manifest must retain its recorded blinding limitation")

    manifest = json.loads((data / required[2]).read_text(encoding="utf-8"))
    final = json.loads((data / required[3]).read_text(encoding="utf-8"))
    values = manifest["results"]
    score = final["data_score"]
    same("LRGxELG amplitude", values["wake_amplitude"], score["amplitude"], problems)
    same("LRGxELG uncertainty", values["hartlap_standard_error"], score["sigma_hartlap"], problems)
    same("LRGxELG nominal z", values["nominal_signed_z"], score["z_signed"], problems)
    same("LRGxELG delta chi2", values["nominal_delta_chi2"], score["delta_chi2_hartlap"], problems)
    for name, final_name in (("empirical_absolute_z", "empirical_abs_z"),
                             ("empirical_sellentin_heavens", "empirical_sellentin_heavens_lr")):
        a, b = values[name], final[final_name]
        if a["n_ge_data"] != b["n_ge_data"]:
            problems.append(name + ": mock tail counts differ")
        same(name + " p", a["p_plus_one"], b["p_plus_one"], problems)
        same(name + " p from counts", b["p_plus_one"],
             (b["n_ge_data"] + 1) / (final["n_mocks"] + 1), problems)
    if manifest["covariance"]["realizations"] != final["n_mocks"]:
        problems.append("LRGxELG final mock count differs from manifest")

    ez = data / "phase7_ezmock_local"
    meta = json.loads((ez / "aggregate_validation_summary.json").read_text(encoding="utf-8"))
    for key, name in (("summary", "summary_ezmock_placebo_covariance.json"),
                      ("oas_covariance", "ezmock_placebo_covariance_oas.csv"),
                      ("sample_covariance", "ezmock_placebo_covariance_sample.csv"),
                      ("vectors", "ezmock_placebo_vectors.csv"),
                      ("window_templates", "ezmock_placebo_window_templates.csv")):
        target = ez / "aggregate" / name
        if not target.is_file() or hashlib.sha256(target.read_bytes()).hexdigest() != meta["aggregate_sha256"][key]:
            problems.append(f"EZmock checksum does not match: {name}")


def main() -> int:
    problems: list[str] = []
    files = tracked_files()
    py, js, sh = validate_sources(files, problems)
    refs = validate_references(files, problems)
    digests = validate_checksums(files, problems)
    validate_reported_results(problems)
    print(f"Checked {py} Python files, {js} JSON files, {sh} shell scripts, "
          f"{refs} local references and {digests} published SHA256 entries.")
    if problems:
        for problem in problems:
            print("ERROR:", problem, file=sys.stderr)
        return 1
    print("Repository integrity checks passed. Survey and physical calculations were not rerun.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
