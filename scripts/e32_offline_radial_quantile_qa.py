#!/usr/bin/env python3
"""E32 synthetic radial-quantile nonidentifiability and analytic sinc enclosure."""
import hashlib
import json
import math
import sys
from fractions import Fraction as Q
from pathlib import Path

HERE=Path(__file__).resolve().parent
OUT=HERE/'e32_offline_radial_quantile_result.json'
PERCENT=[10,25,33,50,67,75,90,95,98,100]
RADII=[Q(i,10) for i in range(1,11)]  # synthetic, dimensionless r100=1

def squantile(shells,p):
    cumulative=Q(0)
    for r,m in sorted(shells):
        cumulative+=m
        if cumulative >= p:
            return r
    raise ValueError('quantile outside normalization')

def sinc(x):
    if x==0: return 1.0
    return math.sin(x)/x

def fourier(shells, kr100):
    return sum(float(m)*sinc(kr100*float(r)) for r,m in shells)

def execute():
    checks=[]
    def verify(name,cond):
        if not cond:raise AssertionError(name)
        checks.append(name)
    a=[];b=[]; last_p=Q(0);last_r=Q(0)
    for pct,r in zip(PERCENT,RADII):
        p=Q(pct,100);w=p-last_p
        verify('positive_bin_mass_'+str(pct),w>0 and r>last_r)
        a.append((r,w))
        b.extend([((last_r+r)/2,w/2),(r,w/2)])
        last_p,last_r=p,r
    verify('both_normalized',sum(m for r,m in a)==sum(m for r,m in b)==1)
    quantiles=[]
    for pct,rad in zip(PERCENT,RADII):
        q=Q(pct,100)
        ra,rb=squantile(a,q),squantile(b,q)
        verify(f'quantile_{pct}_identical',ra==rb==rad)
        quantiles.append({'percent':pct,'radius_over_r100':str(rad)})
    # Negative control: redistribute mass within first bin so that Q(10%) changes.
    bad=list(b)
    first_mid, first_anchor=bad[0],bad[1]
    verify('first_mass_toy_as_expected',first_mid[1]==first_anchor[1]==Q(1,20))
    bad[0]=(first_mid[0],first_mid[1]+first_anchor[1]);bad.pop(1)
    verify('tampered_quantile_rejected',squantile(bad,Q(1,10))!=RADII[0])
    verify('tampered_distribution_still_normalized',sum(m for r,m in bad)==1)
    grid=[]
    for q in (0.1,1.0,3.0,8.0,20.0):
        u_a=fourier(a,q);u_b=fourier(b,q)
        # Quantile information alone bounds the radial characteristic function:
        #  sinc(x)=int_0^1 cos(tx) dt -> |sinc'(x)| <= int_0^1 t dt = 1/2.
        # Each quantile bin has mass w and r in [r_prev,r_i], hence
        # |u - sum w*sinc(q*r_mid)| <= sum w*q*(r_i-r_prev)/4.
        centre=0.;err=0.;last_p=Q(0);last_r=Q(0)
        for pct,r in zip(PERCENT,RADII):
            p=Q(pct,100);w=float(p-last_p)
            rmid=float((last_r+r)/2)
            centre+=w*sinc(q*rmid)
            err+=w*abs(q)*float(r-last_r)/4
            last_p,last_r=p,r
        lo=max(-1.0,centre-err);hi=min(1.0,centre+err)
        verify(f'Fourier_A_inside_enclosure_q_{q}',lo-1e-12<=u_a<=hi+1e-12)
        verify(f'Fourier_B_inside_enclosure_q_{q}',lo-1e-12<=u_b<=hi+1e-12)
        verify(f'Fourier_normalized_modulus_bound_q_{q}',abs(u_a)<=1+1e-12 and abs(u_b)<=1+1e-12)
        grid.append({'k_times_r100':q,'u_shell_A':u_a,'u_shell_B':u_b,
           'absolute_difference':abs(u_b-u_a),'quantile_only_lower':lo,
           'quantile_only_upper':hi,'Lipschitz_radius':err})
    verify('two_profiles_differ_at_positive_k',all(x['absolute_difference']>1e-7 for x in grid))
    verify('toy_high_k_not_a_physical_UV_certification',grid[-1]['k_times_r100']==20)
    report={
       'classification':'E32_SYNTHETIC_MATHEMATICAL_L2_RADIAL_PERCENTILE_NONIDENTIFIABILITY',
       'status':'PASS',
       'no_real_halo_columns_ASDF_FITS_or_observed_data_read':True,
       'units':'dimensionless toy r100=1, q=k*r100; not Abacus halo radii',
       'percentiles':quantiles,
       'mass_in_quantile_bins':[str(Q(p-(PERCENT[i-1] if i else 0),100)) for i,p in enumerate(PERCENT)],
       'A_shell_count':len(a),'B_shell_count':len(b),
       'QA':{'assertions':len(checks),'checks':checks,
           'intentional_mass_reassignment_rejected':True},
       'mathematical_result':{
         'u_spherical':'u(k)=integral sinc(k*r) dM(r)/M',
         'quantile_enclosure':'Let w_i=p_i-p_(i-1), Q(p_i)=r_i, r_0=0. m_i=(r_(i-1)+r_i)/2. Then |u(k)-sum_i w_i sinc(k m_i)| <= |k|/4 sum_i w_i (r_i-r_(i-1)). Clip interval to [-1,1].',
         'derivative_bound':'sinc(x)=int_0^1 cos(tx)dt, |sinc_prime(x)|<=1/2',
         'scope':'radially averaged Fourier mass profile for positive mass relative to a documented L2 centre; does not constrain anisotropic u(vector k), M200c, environment, progenitor, or nonlinear Born correction'},
       'grid':grid,
       'interpretation_stop':{'physical_E28_u(k)':'NOT_IDENTIFIED_FROM_L2_PERCENTILES',
          'real_20_pair_L2_quantile_values':'NOT_READ',
          'UV_total_drag':'UNCERTIFIED',
          'A03':'PHYSICAL_UNCERTIFIED','A04':'BLOCKED','observed_odd':'SEALED'}
    }
    if OUT.exists():
        assert json.loads(OUT.read_text(encoding='utf-8'))==report,'immutable_report_changed'
    else:
        OUT.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('E32_OFFLINE_QA_PASS',len(checks),'assertions')
    for x in grid:
        print('q=',x['k_times_r100'],'u_A=',round(x['u_shell_A'],9),
              'u_B=',round(x['u_shell_B'],9),'abs_delta=',round(x['absolute_difference'],9),
              'quantile_bound=',[round(x['quantile_only_lower'],5),round(x['quantile_only_upper'],5)])
    print('OUTPUT',OUT,'SHA256',hashlib.sha256(OUT.read_bytes()).hexdigest())

if __name__=='__main__':
    try:execute()
    except Exception as e:
        print('E32_OFFLINE_QA_FAIL',type(e).__name__,str(e),flush=True)
        sys.exit(1)
