#!/usr/bin/env python3
"""Posthoc E29 20-pair spatial-selection diagnostic; reads only small JSONs.

No ASDF, FITS, galaxy rows, PID arrays, web access, new selection, or inference.
This is a descriptive QA audit of the first-20-from-256 convenience sample.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import median

BOX_HMPC = 2000.0
EPS = 1e-10


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def min_image(a, b):
    return tuple((float(x)-float(y)) - math.floor(float(x)-float(y)+0.5)
                 for x, y in zip(a, b))


def periodic_cover(values):
    """Shortest circular arc containing all 1-periodic scalar coordinates."""
    s = sorted((float(x)+0.5) % 1.0 for x in values)
    gaps = [(s[(i+1) % len(s)] + (1.0 if i == len(s)-1 else 0.0)) - v
            for i, v in enumerate(s)]
    index = max(range(len(s)), key=lambda i: gaps[i])
    start = s[(index+1) % len(s)]
    width = 1.0 - gaps[index]
    within = [(v-start) % 1.0 for v in s]
    assert abs(max(within)-width) < EPS
    assert all(0.0 <= x < 1.0+EPS for x in within)
    return {
        "minimal_cover_width_unit_box": width,
        "minimal_cover_width_comoving_hinv_mpc": width*BOX_HMPC,
        "fraction_of_full_box_axis": width,
        "minimum_arc_origin_mod_unit_box": start,
        "raw_unit_box_min": min(values),
        "raw_unit_box_max": max(values),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalog', type=Path, default=Path('/mnt/data/e29_20_halo_two_epoch_pair_catalog.json'))
    ap.add_argument('--selection', type=Path, default=Path('/mnt/data/e29_z0950_20_halo_candidates.json'))
    ap.add_argument('--parent-audit', type=Path, default=Path('/mnt/data/e29_20_pair_posthoc_offline_kinematic_audit.json'))
    ap.add_argument('--out', type=Path, default=Path('/mnt/data/e29_20_pair_posthoc_spatial_selection_audit.json'))
    args = ap.parse_args()
    catalog = json.loads(args.catalog.read_text(encoding='utf-8'))
    selection = json.loads(args.selection.read_text(encoding='utf-8'))
    prior = json.loads(args.parent_audit.read_text(encoding='utf-8'))
    pairs = catalog['pairs']
    chosen = selection['selected_halos']
    assert len(pairs) == len(chosen) == 20
    assert catalog['summary']['global_progenitor_certified'] == 0
    assert selection['row_window'] == [0, 256]
    assert selection['selection'] == {'min_N': 100, 'min_npoutA': 5, 'max_halos': 20}
    assert prior['summary']['all_source_mappings_and_arithmetic_checks_pass'] is True
    assert prior['source_sha256'][args.catalog.name] == sha256(args.catalog)
    assert len({p['early_row'] for p in pairs}) == 20
    assert len({p['late_halo_id'] for p in pairs}) == 20
    assert all(not p['global_progenitor_certified'] for p in pairs)
    assert [(p['late_row'],p['late_halo_id'],p['late_N']) for p in pairs] == [
        (q['row'],q['halo_id'],q['N']) for q in chosen]

    displacements = []
    locations = []
    previous = {p['late_halo_id']:p for p in prior['pairs']}
    for p in pairs:
        loc = tuple(float(x) for x in p['late_x_com'])
        early = tuple(float(x) for x in p['early_x_com'])
        disp = min_image(loc, early)
        assert norm(tuple(disp[i]-p['period1_delta'][i] for i in range(3))) < EPS
        assert abs(norm(disp)-p['period1_distance']) < EPS
        assert abs(norm(disp)*BOX_HMPC - previous[p['late_halo_id']]['displacement_comoving_hinv_mpc']) < EPS
        # Periodic-coordinate invariance for each axis independently.
        for i in range(3):
            shifted = list(loc)
            shifted[i] += 1.0
            assert norm(tuple(a-b for a,b in zip(min_image(shifted,early),disp))) < EPS
            shifted = list(early)
            shifted[i] -= 1.0
            assert norm(tuple(a-b for a,b in zip(min_image(loc,shifted),disp))) < EPS
        displacements.append(tuple(x*BOX_HMPC for x in disp))
        locations.append(loc)

    mean = [sum(d[i] for d in displacements)/20 for i in range(3)]
    residual = [tuple(d[i]-mean[i] for i in range(3)) for d in displacements]
    rms = lambda seq: math.sqrt(sum(norm(d)**2 for d in seq)/len(seq))
    raw_width_y = max(p[1] for p in locations)-min(p[1] for p in locations)
    spatial = {ax:periodic_cover([p[i] for p in locations])
               for i,ax in enumerate('xyz')}
    report = {
        'classification':'E29_POSTHOC_OFFLINE_SPATIAL_SELECTION_AND_CENTRE_DISPLACEMENT_QA',
        'source_sha256': {p.name:sha256(p) for p in [args.catalog,args.selection,args.parent_audit]},
        'inputs_no_ASDF_FITS_or_observed_galaxy_rows_read': True,
        'simulation':catalog['simulation'],
        'epochs':{'early':catalog['early_epoch'],'late':catalog['late_epoch']},
        'sample_contract':{
            'late_superslab':0,
            'original_first_rows_window':selection['row_window'],
            'selection':selection['selection'],
            'pair_count':20,
            'early_match_scope_counts':catalog['summary']['early_first_100k_scope'],
            'targeted_spatial_scope_counts':catalog['summary']['targeted_spatial_scope'],
            'globally_certified_progenitor_count':0,
            'not_random_not_representative_population_sample':True,
        },
        'unit_contract':'raw normalized unit-box positions, period 1, BoxSizeHMpc=2000',
        'minimal_periodic_cover_of_late_centres_by_axis':spatial,
        'periodic_wrap_negative_control':{
            'late_y_raw_max_minus_min_unit_box':raw_width_y,
            'late_y_correct_minimal_periodic_arc_unit_box':spatial['y']['minimal_cover_width_unit_box'],
            'integer_period_shift_invariance_all_20_pairs_all_3_axes':'PASS',
        },
        'centre_displacement_descriptive':{
            'mean_vector_comoving_hinv_mpc':mean,
            'norm_of_mean_vector_comoving_hinv_mpc':norm(mean),
            'rms_individual_magnitude_comoving_hinv_mpc':rms(displacements),
            'median_individual_magnitude_comoving_hinv_mpc':median(map(norm,displacements)),
            'rms_after_subtracting_sample_mean_comoving_hinv_mpc':rms(residual),
            'median_after_subtracting_sample_mean_comoving_hinv_mpc':median(map(norm,residual)),
            'norm_sample_mean_over_individual_rms':norm(mean)/rms(displacements),
            'positive_negative_axis_displacement_counts':{
                ax:{'positive':sum(d[i]>0 for d in displacements),'negative':sum(d[i]<0 for d in displacements)}
                for i,ax in enumerate('xyz')
            },
        },
        'interpretation_stop':{
            'selection':'The first 20 eligible halos in first 256 superslab-000 rows form a narrow x/y stripe, not a representative or independent halo population.',
            'displacement':'The sample-average centre displacement is descriptive, not a neutrino-CDM wind measurement or a halo acceleration.',
            'mass_and_profile':'L1 N is not M200c; no same-object u_src(k,z), u_test(k,zobs), environment, or certified UV recoil.',
            'progenitors':'All 20 matches are unique only in their previously tested early-catalog windows, not a globally validated merger tree.',
            'observations':'No observed eBOSS galaxy rows or 24D odd vector used; A03 PHYSICAL_UNCERTIFIED, A04 BLOCKED.',
        },
    }
    tmp = args.out.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    tmp.replace(args.out)
    print('E29_OFFLINE_SPATIAL_SELECTION_AUDIT_PASS')
    print('SELECTION',report['sample_contract'])
    print('WIDTHS_HINV_MPC',*[round(spatial[ax]['minimal_cover_width_comoving_hinv_mpc'],6) for ax in 'xyz'])
    print('MEAN_VECTOR_HINV_MPC',*[round(x,8) for x in mean])
    print('MEAN_NORM_RMS_RESIDUAL',round(norm(mean),8),round(rms(displacements),8),round(rms(residual),8))
    print('Y_PERIODIC_WRAP_NEGATIVE_CONTROL_PASS',round(raw_width_y,8),round(spatial['y']['minimal_cover_width_unit_box'],8))
    print('REPORT',args.out)
    print('REPORT_SHA256',sha256(args.out))

if __name__ == '__main__':
    main()
