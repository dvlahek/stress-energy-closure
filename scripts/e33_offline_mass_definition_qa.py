#!/usr/bin/env python3
"""E33 offline flat matter+Lambda critical-density pseudo-evolution toy."""
import hashlib
import json
import math
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'e29_20_pair_posthoc_offline_kinematic_audit.json'
OUT=HERE/'e33_offline_mass_definition_result.json'

def execute():
    checks=[]
    def verify(name,cond):
        if not cond:raise AssertionError(name)
        checks.append(name)
    parent=json.loads(SOURCE.read_text(encoding='utf-8'))
    verify('parent_is_20_pair_previous_offline_audit',parent['summary']['pairs']==20 and parent['inputs_no_ASDF_FITS_or_observed_galaxy_rows_read'])
    p=parent['cosmology_time_approximation']
    ze=p['redshift_early_actual_header'];zl=p['redshift_late_actual_header']
    om=p['Omega_m'];ol=p['Omega_Lambda']
    verify('correct_time_order',ze>zl>0)
    verify('flat_background_approximation',abs(om+ol-1)<1e-12 and 0<om<1)
    def e2(z,mo=om,la=ol):return mo*(1+z)**3+la
    ec=e2(ze);lc=e2(zl);rat=lc/ec
    verify('critical_density_decreases',0<rat<1)
    verify('same_redshift_no_change',e2(zl)/e2(zl)==1)
    verify('de_sitter_null',abs(e2(zl,0,1)/e2(ze,0,1)-1)<1e-15)
    examples=[]
    for gamma in (2.,2.5):
        rr=rat**(-1/gamma)
        mr=rr**(3-gamma)
        verify(f'static_profile_r200_moves_{gamma}',rr>1)
        verify(f'static_profile_M200_increases_{gamma}',mr>1)
        verify(f'analytic_mass_ratio_{gamma}',abs(mr-rat**(-(3-gamma)/gamma))<1e-13)
        examples.append({'gamma':gamma,'physical_profile':'static_power_law_toy',
           'R200c_late_over_early':rr,'M200c_late_over_early':mr,
           'apparent_M200c_growth_percent':(mr-1)*100,
           'physical_accretion_percent_by_construction':0.0})
    report={'classification':'E33_POSTHOC_OFFLINE_M200C_REFERENCE_DENSITY_PSEUDOEVOLUTION_TOY',
       'input_file':SOURCE.name,'input_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
       'no_asdf_fits_observed_data_or_new_halo_column_read':True,
       'source_epoch':{'z_early':ze,'z_late':zl,'approx_H0_km_s_Mpc':p['H0_km_s_Mpc'],
          'Omega_m0':om,'Omega_Lambda0':ol,'cosmology_approximation':p['assumption']},
       'E2_early':ec,'E2_late':lc,'rho_crit_late_over_early':rat,
       'rho_crit_percent_change':(rat-1)*100,
       'example_static_halos':examples,
       'QA':{'status':'PASS','assertions':len(checks),'checks':checks},
       'scope':'Purely hypothetical fixed physical density rho(r) proportional r^-gamma; do not transfer example masses to measured L1 N or actual M200c',
       'physical_stop':{'L1_N_as_M200c':'NOT_PERMITTED',
         'physical_accretion_from_two_snapshot_N':'NOT_IDENTIFIED',
         'M200c_from_real_20_halo_profiles':'NOT_MEASURED',
         'A03':'PHYSICAL_UNCERTIFIED','A04':'BLOCKED','observed_odd':'SEALED'}}
    if OUT.exists():
        assert json.loads(OUT.read_text(encoding='utf-8'))==report,'immutable_report_changed'
    else:
        OUT.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('E33_OFFLINE_QA_PASS',len(checks),'assertions')
    print('RHOCRIT_LATE_OVER_EARLY',format(rat,'.12f'),'CHANGE_PERCENT',format((rat-1)*100,'.6f'))
    for ex in examples:
        print('STATIC_GAMMA',ex['gamma'],'APPARENT_M200C_GROWTH_PERCENT',format(ex['apparent_M200c_growth_percent'],'.6f'))
    print('OUTPUT',OUT,'SHA256',hashlib.sha256(OUT.read_bytes()).hexdigest())

if __name__=='__main__':
    try:execute()
    except Exception as e:
        print('E33_OFFLINE_QA_FAIL',type(e).__name__,str(e),flush=True)
        sys.exit(1)
