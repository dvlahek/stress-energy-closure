#!/usr/bin/env python3
"""Offline audit of user-supplied E29 JSON only; opens no ASDF or FITS."""
import hashlib
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FN = {
    'catalog': 'e29_20_halo_two_epoch_pair_catalog.json',
    'selection': 'e29_z0950_20_halo_candidates.json',
    'window': 'e29_20_halo_pid_extended_100k.json',
    'periodic': 'e29_20_halo_spatial_audit_period1.json',
    'spatial': 'e29_three_unresolved_spatial_candidates.json',
    'targeted': 'e29_three_unresolved_targeted_pid_matches.json',
}
raw = {k: (ROOT / f).read_bytes() for k, f in FN.items()}
src = {k: json.loads(v) for k, v in raw.items()}
cat, selection, window, periodic, spatial, targeted = [src[k] for k in FN]
pairs = cat['pairs']
assert len(pairs) == len({p['late_halo_id'] for p in pairs}) == 20
assert len({p['early_row'] for p in pairs}) == 20
assert cat['summary']['selected_late_halos'] == 20
assert spatial['complete'] and spatial['rows_scanned'] == 11608840
assert len(spatial['targets']) == len(targeted['targets']) == 3
win_by_id = {m['late_halo_id']: m for m in window['matches']}
per_by_id = {m['late_halo_id']: m for m in periodic['halos']}
tgt_by_id = {m['late_halo_id']: m for m in targeted['targets']}
selected_by_id = {m['halo_id']: m for m in selection['selected_halos']}
assert len(selected_by_id) == len(pairs) == 20

for p in pairs:
    late_id = p['late_halo_id']
    initial = selected_by_id[late_id]
    assert (p['late_row'],p['late_N'],p['late_pid_count']) == (initial['row'],initial['N'],initial['npoutA'])
    if p['match_scope'] == 'EARLY_FIRST_100000_ROWS_N_GE_50':
        source = win_by_id[late_id]
        assert source['qualifying_candidate_count'] == 1
        assert source['status'] == 'SINGLE_QUALIFYING_CANDIDATE_IN_100K_WINDOW'
        candidate = source['candidates'][0]
        assert per_by_id[late_id]['early_halo_id'] == p['early_halo_id']
        assert per_by_id[late_id]['early_N'] == p['early_N']
    else:
        assert p['match_scope'] == 'SPATIAL_RADIUS_0.001_ALL_EARLY_SUPERSLAB_000_N_GE_50'
        source = tgt_by_id[late_id]
        assert source['qualifying_candidates_ge5_shared'] == 1
        assert source['status'] == 'SINGLE_IN_SPATIAL_WINDOW'
        candidate = next(c for c in source['candidates'] if c['shared_pid_count'] >= 5)
        assert candidate['early_halo_id'] == p['early_halo_id']
    for key in ('early_halo_id', 'early_row', 'early_N', 'shared_pid_count'):
        if key in candidate:
            assert candidate[key] == p[key], (late_id, key)
    assert p['shared_pid_count'] >= 5
    assert p['global_progenitor_certified'] is False
    assert p['late_N'] >= 100 and p['early_N'] >= 50
    delta = [x-y for x,y in zip(p['late_x_com'],p['early_x_com'])]
    wrap = [d-math.floor(d+0.5) for d in delta]
    assert max(abs(x-y) for x,y in zip(wrap,p['period1_delta'])) < 1.e-12
    assert abs(math.dist((0,0,0),wrap)-p['period1_distance']) < 1.e-12
    assert p['delta_N'] == p['late_N']-p['early_N']
    assert abs(p['relative_delta_N']-p['delta_N']/p['early_N']) < 1.e-12
    assert abs(p['late_pid_overlap']-p['shared_pid_count']/p['late_pid_count'])<1.e-12

# From uploaded Abacus YAML-only header; no new ASDF array/header read.
z_early=1.027997062494471
z_late=0.952838237036305
h=0.6736
H0=67.36
omega_m=0.315192
omega_l=0.684808
box_hmpc=2000.0
mpc_km=3.0856775814913673e19
sec_gyr=365.25*86400*1e9

def simpson(f,a,b,n=4096):
    assert n%2==0
    step=(b-a)/n
    return step/3*(f(a)+f(b)+4*sum(f(a+step*i) for i in range(1,n,2))+2*sum(f(a+step*i) for i in range(2,n,2)))

a_early=1/(1+z_early);a_late=1/(1+z_late)
H_inv_seconds=mpc_km/H0
delta_t_gyr = H_inv_seconds/sec_gyr*simpson(lambda a: 1/(a*math.sqrt(omega_m/a**3+omega_l)), a_early,a_late)
# Proper-velocity proxy computed at mean scale factor, not a measured halo peculiar velocity.
a_mid=(a_early+a_late)/2
pair_rows=[]
for p in pairs:
    comov_hmpc=p['period1_distance']*box_hmpc
    comov_mpc=comov_hmpc/h
    speed_kms=a_mid*comov_mpc*mpc_km/(delta_t_gyr*sec_gyr)
    pair_rows.append({
        'late_halo_id':p['late_halo_id'],
        'early_halo_id':p['early_halo_id'],
        'shared_pid_count':p['shared_pid_count'],
        'late_pid_overlap':p['late_pid_overlap'],
        'early_N':p['early_N'],
        'late_N':p['late_N'],
        'relative_delta_N':p['relative_delta_N'],
        'displacement_comoving_hinv_mpc':comov_hmpc,
        'displacement_proper_at_mean_a_Mpc':a_mid*comov_mpc,
        'finite_difference_centre_speed_proxy_km_s':speed_kms,
        'scope':p['match_scope'],
    })

report={
 'classification':'E29_POSTHOC_OFFLINE_20_LOCAL_PID_PAIR_CONSISTENCY_AND_KINEMATIC_PROXY',
 'source_sha256':{FN[k]:hashlib.sha256(v).hexdigest() for k,v in raw.items()},
 'inputs_no_ASDF_FITS_or_observed_galaxy_rows_read':True,
 'status':'POSTHOC_DESCRIPTIVE_ONLY_NOT_GLOBAL_PROGENITOR_OR_M200C_HISTORY',
 'cosmology_time_approximation':{
   'redshift_early_actual_header':z_early,'redshift_late_actual_header':z_late,
   'H0_km_s_Mpc':H0,'Omega_m':omega_m,'Omega_Lambda':omega_l,
   'assumption':'flat matter+Lambda approximation to snapshot spacing; not exact CLASS/Abacus time integral with massive neutrinos',
   'time_spacing_Gyr_approx':delta_t_gyr,
   'a_mean':a_mid,
 },
 'source_position_units':'raw unit box (official AbacusSummit data-products specification); BoxSizeHMpc=2000; minimal-image period=1',
 'summary':{
   'pairs':20,
   'window_100k_pairs':sum('FIRST_100000' in p['match_scope'] for p in pairs),
   'spatial_radius_000_pairs':sum('SPATIAL_RADIUS' in p['match_scope'] for p in pairs),
   'all_source_mappings_and_arithmetic_checks_pass':True,
   'global_progenitor_certified':0,
   'min_shared_PIDs':min(p['shared_pid_count'] for p in pairs),
   'median_late_PID_overlap':statistics.median(p['late_pid_overlap'] for p in pairs),
   'median_relative_N_change':statistics.median(p['relative_delta_N'] for p in pairs),
   'displacement_comoving_hinv_mpc_min':min(p['displacement_comoving_hinv_mpc'] for p in pair_rows),
   'displacement_comoving_hinv_mpc_median':statistics.median(p['displacement_comoving_hinv_mpc'] for p in pair_rows),
   'displacement_comoving_hinv_mpc_max':max(p['displacement_comoving_hinv_mpc'] for p in pair_rows),
   'finite_difference_centre_speed_proxy_km_s_median':statistics.median(p['finite_difference_centre_speed_proxy_km_s'] for p in pair_rows),
   'finite_difference_centre_speed_proxy_km_s_max':max(p['finite_difference_centre_speed_proxy_km_s'] for p in pair_rows),
 },
 'interpretation_stop':{
   'selection':'20 first qualifying late halos from first 256 rows of superslab 000; not a representative population sample',
   'N':'CompaSO L1 assigned particle count, not M200c or physical accretion history',
   'PID':'3 percent subsample A; local unique within tested window not global merger-tree confirmation',
   'speed':'snapshot-centre finite difference, not neutrino-CDM relative wind, instantaneous peculiar velocity, or E28 v_h=+200 km/s validation',
   'profile':'no same-object normalized physical source/test u(k,z) or certified UV force',
   'A03_A04':'A03 PHYSICAL_UNCERTIFIED; A04 BLOCKED; observed 24D eBOSS odd SEALED'
 },
 'pairs':pair_rows,
}
out=ROOT/'e29_20_pair_posthoc_offline_kinematic_audit.json'
out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print('OFFLINE_E29_20_PAIR_AUDIT_PASS')
print('NO_ABACUS_ASDF_OR_FITS_OPENED')
print('SOURCE_MAPS_RECHECKED',len(pair_rows))
print('DELTA_TIME_GYR_APPROX',round(delta_t_gyr,7))
print('DISPLACEMENT_HINV_MPC_MIN_MED_MAX',*[round(report['summary'][f'displacement_comoving_hinv_mpc_{k}'],6) for k in ('min','median','max')])
print('MEDIAN_CENTRE_SPEED_PROXY_KM_S',round(report['summary']['finite_difference_centre_speed_proxy_km_s_median'],1))
print('REPORT',out, 'BYTES',out.stat().st_size)
