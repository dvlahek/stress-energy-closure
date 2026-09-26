#!/usr/bin/env python3
"""Reproduce SHA-pinned realistic EZmock0001 cross pairs, then project mock xi_l.

Only fixed 0001 mock galaxies and their same-realization tracer/cap randoms.
All eight complete compressed SHA256s and all sixteen original forward/reverse
histogram SHAs must match the ARCHIVED first pilot before any new projection.
The newly reported ell=1,3 are descriptive SINGLE-MOCK values, not the
observed odd-sector statistic, a calibrated physical null or a covariance.
"""
from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np

from audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot import (
    ROOT, PILOT as PARENT_PROTOCOL, ORDER, ORDER_KEYS,
    source_manifest_gate, verified_eight_sources, sample_catalogue,
    oriented_pair_terms,
)
from check_eboss_cross_ls_synthetic import cross_landy_szalay, LABELS
from eboss_dr16_fiducial import PRIMARY_GEOMETRY, comoving_mpc_over_h
from audit_eboss_dr16_rr_pair_closure import mirrored_closure

PROTOCOL = ROOT / "source_data/eboss_dr16_ezmock0001_galaxy_odd_projection_protocol_2026-09-26.json"
REPORT = ROOT / "source_data/eboss_dr16_ezmock0001_galaxy_cross_ls_pilot_report_2026-09-26.json"
MANIFEST = ROOT / "source_data/eboss_dr16_ezmock0001_galaxy_cross_ls_pilot_uploaded_manifest_2026-09-26.json"
PASSED = "EZMOCK0001_MOCK_GALAXY_ODD_PROJECTION_DESCRIPTIVE_ONLY"
FAILED = "EZMOCK0001_MOCK_GALAXY_ODD_PROJECTION_INCOMPLETE_STOP"
S = np.asarray([20., 40., 60., 80., 100., 120., 140.])
MU = np.linspace(-1 - 1e-7, 1 + 1e-7, 25)


def git_blob_sha(raw):
    return hashlib.sha1(
        ("blob " + str(len(raw))).encode("ascii") + bytes([0]) + raw
    ).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def check_protocol(p):
    fixed = read_json(MANIFEST)
    prior = read_json(PARENT_PROTOCOL)
    raw = REPORT.read_bytes()
    local = (ROOT / p["parent_pilot_local"]).read_bytes()
    if (
        raw != local
        or len(raw) != fixed["exact_uploaded_bytes"]
        or hashlib.sha256(raw).hexdigest() != fixed["exact_uploaded_sha256"]
        or hashlib.sha256(raw).hexdigest() != p["parent_pilot_report_exact_sha256"]
        or git_blob_sha(raw) != fixed["exact_uploaded_git_blob_sha1"]
        or git_blob_sha(raw) != p["parent_pilot_report_git_blob_sha1"]
        or p["parent_pilot_report_archive"] != str(REPORT.relative_to(ROOT))
        or p["parent_pilot_manifest"] != str(MANIFEST.relative_to(ROOT))
        or p["mock_id"] != prior["mock_id"]
        or p["mock_id"] != 1
        or p["caps"] != prior["caps"]
        or p["tracers"] != prior["tracers"]
        or p["sample_count_per_tracer_cap"] != {
            "galaxy_D":prior["pilot_counts_per_cap_tracer"]["mock_galaxy_D"],
            "matched_random_R":prior["pilot_counts_per_cap_tracer"]["mock_random_R"],
        }
        or p["redshift"] != prior["exact_fixed_highz_bin"]
        or p["s_edges_mpc_h"] != S.tolist()
        or prior["geometry"]["s_edges_mpc_h"] != S.tolist()
        or p["mu_edges"] != "np.linspace(-1 - 1e-7, 1 + 1e-7, 25)"
        or p["theta_min_deg"] != prior["geometry"]["theta_min_deg"]
        or p["fiducial"] != PRIMARY_GEOMETRY
        or p["legendre_orders"] != [0,1,2,3]
        or p["observed_galaxy_rows_read"] is not False
        or p["observed_odd_data_vector_read"] is not False
        or p["new_science_selection_applied"] is not False
        or p["mock_covariance_computed"] is not False
    ):
        raise ValueError("Frozen parent result, fixed geometry or odd projection protocol changed")
    old = json.loads(raw)
    if (
        old["status"] != fixed["report_status"]
        or old["pilot_mock_id"] != 1
        or old["errors"] != []
        or old["mock_galaxy_rows_read"] is not True
        or old["mock_random_rows_read"] is not True
        or old["observed_galaxy_rows_read"] is not False
        or old["observed_odd_data_vector_read"] is not False
        or old["new_science_selection_applied"] is not False
        or old["mock_covariance_computed"] is not False
        or old["exact_uploaded_raw_report_sha256"] != fixed["parent_source_sha_report_sha256"]
        or old["source_galaxy_full_compressed_SHA256"] != fixed["source_galaxy_sha256_4"]
        or old["matched_source_random_full_compressed_SHA256"] != fixed["matched_random_sha256_4"]
        or [case["cap"] for case in old["cases"]] != ["NGC","SGC"]
    ):
        raise ValueError("The uploaded original matched mock0001 pair report fails pre-projection gate")
    for case in old["cases"]:
        if (
            case["cross_ls"]["RR_supported_s_mu_cells"] != 144
            or case["cross_ls"]["RR_total_s_mu_cells"] != 144
            or set(case["cross_ls"]["forward_pair_terms"]) != set(LABELS)
            or set(case["cross_ls"]["reverse_pair_terms"]) != set(LABELS)
        ):
            raise ValueError("Original fixed mock0001 full RR-support invariant changed")
    return old, prior


def legendre_integral_at(x, ell):
    if ell == 0:
        return x
    if ell == 1:
        return 0.5*x*x
    if ell == 2:
        return 0.5*(x*x*x-x)
    if ell == 3:
        return (5./8.)*x**4 - (3./4.)*x*x
    raise ValueError("Only ell=0,1,2,3 predeclared")


def integration_weights(muedges, ell):
    low = np.maximum(-1., muedges[:-1])
    high = np.minimum(1., muedges[1:])
    if not np.all(high > low):
        raise ValueError("Frozen signed-mu intervals do not meet physical [-1,1]")
    weights = (2*ell+1)/2. * (
        legendre_integral_at(high,ell) - legendre_integral_at(low,ell)
    )
    if not np.allclose(weights, ((-1)**ell)*weights[::-1],
                       rtol=0, atol=1e-13):
        raise ValueError("Physical-domain mu integration breaks Legendre parity")
    return np.asarray(weights,dtype="f8")


def project(xi,support,muedges,ells):
    if (
        xi.shape != (6,24) or support.shape != xi.shape
        or not np.all(support)
        or not np.isfinite(xi).all()
        or not np.all(np.diff(muedges)>0)
    ):
        raise ValueError("All 144 signed-mu cells must be positive RR-supported, with finite xi")
    return {
        str(ell): {
            "values_by_fixed_s_bin": (xi @ integration_weights(muedges,ell)).tolist(),
            "integration_weight_array_sha256": hashlib.sha256(
                np.ascontiguousarray(integration_weights(muedges,ell)).tobytes()
            ).hexdigest(),
        }
        for ell in ells
    }


def verify_pair_reproduction(old, cap, label, h, info, orient):
    item=next(v for v in old["cases"] if v["cap"] == cap)
    expected = item["cross_ls"][orient+"_pair_terms"][label]
    got=hashlib.sha256(np.ascontiguousarray(h).tobytes()).hexdigest()
    if (
        got != expected["weighted_histogram_SHA256"]
        or info["accepted_pairs"] != expected["accepted_pairs"]
        or not np.isclose(info["pair_normalization"],
                          expected["independently_normalized_pair_weight"],
                          rtol=1e-13,atol=0)
    ):
        raise ValueError("Archived mock0001 weighted histogram is NOT reproduced exactly: "
                         + cap + "/" + orient + "/" + label)


def audit(p, *, no_download):
    old, previous_protocol = check_protocol(p)
    orig, source_report, source_manifest, random_manifest, highz = (
        source_manifest_gate(previous_protocol)
    )
    # This gate rehashes all FOUR complete galaxy gzip bytes first, then
    # all FOUR complete companion random gzip bytes. No FITS rows beforehand.
    paths = verified_eight_sources(
        orig,source_report,source_manifest,random_manifest,
        no_download=no_download
    )
    output={
        "status":PASSED,
        "original_pilot_report_exact_sha256":p["parent_pilot_report_exact_sha256"],
        "original_sample_counts":p["sample_count_per_tracer_cap"],
        "fixed_s_edges_mpc_h":p["s_edges_mpc_h"],
        "fixed_signed_mu_edge_definition":p["mu_edges"],
        "legendre_orders":p["legendre_orders"],
        "source_compressed_SHA_and_prior_pair_histograms_reverified":True,
        "caps":{},
        "mock_galaxy_and_matched_random_rows_read":True,
        "observed_galaxy_rows_read":False,
        "observed_odd_data_vector_read":False,
        "new_science_selection_applied":False,
        "not_a_physical_odd_null_or_detection":True,
        "not_a_mock_covariance":True,
        "full_physical_LRG_ELG_mask_or_window_certified":False,
        "errors":[],
    }
    distance=lambda z:comoving_mpc_over_h(z, PRIMARY_GEOMETRY)
    mapping={"D1D2":"D1D2","D1R2":"R1D2",
             "R1D2":"D1R2","R1R2":"R1R2"}
    for cap in p["caps"]:
        prior_case=next(x for x in old["cases"] if x["cap"]==cap)
        cats={}
        for tracer in p["tracers"]:
            key=cap+"/"+tracer
            for role in ("dat","ran"):
                cat,metadata=sample_catalogue(
                    paths[key+"/"+role],
                    expected_rows=(
                        source_manifest["per_source_frozen_from_uploaded_report"][key][
                            "expected_mock_galaxy_header_rows_from_prior_header_audit"
                        ] if role=="dat" else int(random_manifest[key]["header_rows"])
                    ),
                    cap=cap,tracer=tracer,role=role,
                    expected_highz=(
                        int(highz["caps"][cap]["mock_randoms"][tracer.removeprefix("eBOSS_")])
                        if role=="ran" else None
                    ),
                    p=previous_protocol
                )
                expected=prior_case["input_sample_diagnostics"][tracer+"_"+role]
                if (
                    metadata["selected_array_SHA256"] != expected["selected_array_SHA256"]
                    or metadata["selected_rows"] != expected["selected_rows"]
                    or metadata["eligible_highz_rows_after_fixed_weight_gate"] !=
                       expected["eligible_highz_rows_after_fixed_weight_gate"]
                ):
                    raise ValueError("Original 0001 sample is not reproducible before odd projection: "
                                     + key+"/"+role)
                cats[tracer,role]=cat
        lrg_d,lrg_r = cats["eBOSS_LRG","dat"],cats["eBOSS_LRG","ran"]
        elg_d,elg_r = cats["eBOSS_ELG","dat"],cats["eBOSS_ELG","ran"]
        fwd,fnorm,fmeta=oriented_pair_terms(lrg_d,elg_d,lrg_r,elg_r,
                                             distance=distance)
        rev,rnorm,rmeta=oriented_pair_terms(elg_d,lrg_d,elg_r,lrg_r,
                                             distance=distance)
        for label in LABELS:
            verify_pair_reproduction(old,cap,label,fwd[label],fmeta[label],"forward")
            verify_pair_reproduction(old,cap,label,rev[label],rmeta[label],"reverse")
            mirror=mirrored_closure(fwd[label],rev[mapping[label]],
                         {"pair_normalization":fnorm[label],
                          "accepted_pairs":fmeta[label]["accepted_pairs"]},
                         {"pair_normalization":rnorm[mapping[label]],
                          "accepted_pairs":rmeta[mapping[label]]["accepted_pairs"]},
                         MU)
            if not mirror["closure_passed"]:
                raise ValueError("Original mock cross LS mirror pair closure no longer valid")
        xi,support=cross_landy_szalay(fwd,fnorm)
        xi_rev,support_rev=cross_landy_szalay(rev,rnorm)
        if (
            not np.array_equal(support,support_rev[:,::-1])
            or not np.all(support) or not np.all(support_rev)
        ):
            raise ValueError("All predeclared signed-mu bins must be RR supported")
        proj=project(xi,support,MU,p["legendre_orders"])
        revproj=project(xi_rev,support_rev,MU,p["legendre_orders"])
        parity={}
        for ell in p["legendre_orders"]:
            a=np.asarray(proj[str(ell)]["values_by_fixed_s_bin"])
            b=np.asarray(revproj[str(ell)]["values_by_fixed_s_bin"])
            residual=float(np.max(np.abs(a-(-1)**ell*b)))
            if not np.isfinite(residual) or residual >= 1e-10:
                raise ValueError("Odd/even projected tracer reversal failed for ell="+str(ell))
            parity[str(ell)]={"forward_reverse_max_abs_parity_residual":residual}
        output["caps"][cap]={
            "RR_supported_s_mu_cells":int(np.count_nonzero(support)),
            "RR_total_s_mu_cells":int(support.size),
            "xi_grid_SHA256":hashlib.sha256(np.ascontiguousarray(xi).tobytes()).hexdigest(),
            "mock_galaxy_xi_multipoles_forward":proj,
            "mock_galaxy_xi_multipoles_reverse":revproj,
            "forward_reverse_projection_parity":parity,
            "single_fixed_mock_and_finite_sample_NOT_inference":True,
        }
        print("MOCK_ODD_PROJECTION_CAP",cap,
              "P1",json.dumps(proj["1"]["values_by_fixed_s_bin"]),
              "P3",json.dumps(proj["3"]["values_by_fixed_s_bin"]),flush=True)
    if set(output["caps"])!=set(p["caps"]):
        raise ValueError("Both fixed mock0001 caps required")
    return output


def self_test(p):
    check_protocol(p)
    xi=np.ones((6,24),dtype="f8")
    proj=project(xi,np.ones((6,24),dtype=bool),MU,[0,1,2,3])
    if not np.allclose(proj["0"]["values_by_fixed_s_bin"],1.0,atol=1e-14):
        raise AssertionError("Constant xi should have unit monopole")
    for ell in (1,2,3):
        if not np.allclose(proj[str(ell)]["values_by_fixed_s_bin"],0.0,atol=1e-12):
            raise AssertionError("Constant synthetic xi has unexpected higher multipole")
    odd=np.broadcast_to(np.linspace(-1,1,24),(6,24)).copy()
    ref=project(odd,np.ones_like(odd,dtype=bool),MU,[1,3])
    reverse=project(odd[:,::-1],np.ones_like(odd,dtype=bool),MU,[1,3])
    if not np.allclose(ref["1"]["values_by_fixed_s_bin"],
                       -np.asarray(reverse["1"]["values_by_fixed_s_bin"]),atol=1e-12):
        raise AssertionError("Synthetic odd projected parity was not preserved")
    support=np.ones((6,24),dtype=bool)
    support[0,0]=False
    try:
        project(xi,support,MU,[1,3])
    except ValueError:
        pass
    else:
        raise AssertionError("Unsupported mock xi cell incorrectly projected")
    print("EBOSS_EZMOCK0001_GALAXY_ODD_PROJECTION_SYNTHETIC_SELF_TEST_OK",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    parser.add_argument("--no-download",action="store_true")
    args=parser.parse_args()
    p=read_json(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0
    outfile=ROOT/p["local_report"]
    try:
        result=audit(p,no_download=args.no_download)
    except Exception as exc:
        result={"status":FAILED,"errors":[str(exc)],
                "mock_galaxy_or_random_rows_may_have_been_read":True,
                "observed_galaxy_rows_read":False,
                "observed_odd_data_vector_read":False,
                "new_science_selection_applied":False,
                "mock_covariance_computed":False}
    outfile.parent.mkdir(parents=True,exist_ok=True)
    temp=outfile.with_suffix(".tmp.json")
    temp.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    temp.replace(outfile)
    print("EZMOCK0001_GALAXY_ODD_PROJECTION",result["status"],flush=True)
    print("REPORT",outfile,flush=True)
    if result["errors"]:
        print("ERRORS",*result["errors"],sep="\n",flush=True)
        return 2
    print("COMPLETED_CAPS",len(result["caps"]),flush=True)
    print("OBSERVED_ODD_DATA_READ",result["observed_odd_data_vector_read"],flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
