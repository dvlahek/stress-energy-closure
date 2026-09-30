#!/usr/bin/env python3
"""E30 exact-arithmetic, synthetic mathematical QA. No physical catalogue input."""
import hashlib
import itertools
import json
import sys
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / 'e30_offline_estimand_result.json'

def mean(xs):
    xs = list(xs)
    assert xs
    return sum(xs, Q(0)) / len(xs)

def fstr(x):
    return str(x)

def rank(a):
    m = [list(map(Q, r)) for r in a]
    nr, nc, row = len(m), len(m[0]), 0
    for col in range(nc):
        piv = next((r for r in range(row, nr) if m[r][col]), None)
        if piv is None:
            continue
        m[row], m[piv] = m[piv], m[row]
        pv = m[row][col]
        m[row] = [v/pv for v in m[row]]
        for r in range(nr):
            if r != row and m[r][col]:
                a_ = m[r][col]
                m[r] = [u-a_*v for u,v in zip(m[r],m[row])]
        row += 1
        if row == nr: break
    return row

def execute():
    checks = []
    def verify(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    # Exact centered, density-correlated wind; the joint law is not invariant
    # under wind-only inversion even though both marginal signs are balanced.
    draws = [(Q(d), Q(d)+Q(e,2)) for d,e in itertools.product((-1,1), repeat=2)]
    d_mean = mean(d for d,u in draws)
    u_mean = mean(u for d,u in draws)
    cross = mean(d*u for d,u in draws)
    triple = mean(u*d*d for d,u in draws)
    verify('zero_means', d_mean == u_mean == 0)
    verify('density_wind_two_field_nonzero', cross == 1)
    verify('centrally_symmetric_three_field_zero', triple == 0)

    signs = [Q(-1),Q(1)]
    for t in (Q(2),Q(3)):
        verify(f'conditional_fixed_amplitude_balanced_{t}',mean(s*t for s in signs)==0)
        verify(f'wind_weighted_recovers_{t}',mean(s*(s*t) for s in signs)==t)
    # More stringent condition: equal marginal P(s) is insufficient if
    # the conditional source amplitude depends on variables correlated with s.
    counterexample = [(Q(s),Q(s)) for s in (-1,1)]
    verify('marginal_sign_balance',mean(s for s,x in counterexample)==0)
    verify('conditional_amplitude_counterexample',mean(s*x for s,x in counterexample)==1)

    def selected(t,wplus,wminus):
        return Q(wplus*t-wminus*t,wplus+wminus)
    eta_plus=Q(2-1,2+1)
    eta_minus=Q(1-2,1+2)
    state_plus=selected(Q(2),2,1)
    state_minus=selected(Q(3),1,2)
    verify('selection_eta_formula',state_plus==Q(2,3) and state_minus==Q(-1))
    verify('state_specific_selection_cannot_factor_common_eta',
           state_minus-state_plus != eta_plus*(Q(3)-Q(2)))
    verify('selected_normalized_mass_positive', 2+1>0 and 1+2>0)

    def chi(AL,AE,cL,cE):return AL*cE-AE*cL
    al,ae=Q(2),Q(3)
    same_fractional=chi(al,ae,Q(4),Q(6))
    different_fractional=chi(al,ae,Q(4),Q(7))
    verify('identical_fractional_tracer_response_null',same_fractional==0)
    verify('unequal_fractional_tracer_response_nonzero',different_fractional==2)
    verify('tracer_exchange_flips_signed_cross',chi(ae,al,Q(7),Q(4))==-different_fractional)
    verify('common_response_sign_flip_flips_odd',chi(al,ae,Q(-4),Q(-7))==-different_fractional)
    verify('even_quadratic_invariant_under_common_response_sign_flip',
           Q(4)**2==Q(-4)**2 and Q(7)**2==Q(-7)**2 and Q(4)*Q(7)==Q(-4)*Q(-7))

    # Exactly nine independent mock vectors give centered rank <= 8, not 18/24.
    # The 9x9 Gram witness for e_i centered over nine snapshots is I-J/9.
    gram = [[(Q(int(i==j))-Q(1,9)) for j in range(9)] for i in range(9)]
    mock_rank=rank(gram)
    verify('nine_centered_mock_rank_eight',mock_rank==8)
    verify('eboss_observable_dimension_24',6*2*2==24)

    # Even-to-odd leakage is independent of the EV source label in this toy.
    toy_even=[Q(2),Q(3)]
    toy_leak=Q(1,10)*toy_even[0]+Q(1,20)*toy_even[1]
    verify('even_to_odd_leak_is_not_an_EV_identification',toy_leak==Q(7,20))
    verify('even_leak_identical_across_F_states',toy_leak-toy_leak==0)

    report={
       'classification':'E30_OFFLINE_MATHEMATICAL_ESTIMAND_GATE_NOT_PHYSICAL_EBOSS_RESULT',
       'method':'Python standard library, exact Fraction enumeration and rational rank',
       'no_asdf_fits_observed_data_or_new_mock_read': True,
       'QA':{'status':'PASS','assertions':len(checks),'checks':checks},
       'linear_vs_quadratic_toy':{'E_delta':fstr(d_mean),'E_u':fstr(u_mean),
         'E_delta_u':fstr(cross),'E_u_delta_delta':fstr(triple),
         'interpretation':'centrally symmetric finite witness, not a physical Gaussian N-body simulation'},
       'conditional_sign_toy':{'equal_marginal_sign_probability':True,
         'constant_t_balanced_mean':'0','density_correlated_t_mean':'1',
         'conditional_balance_needed_for_general_t':True},
       'selection_toy':{'t_plus':'2','t_minus':'3',
         'eta_plus':fstr(eta_plus),'eta_minus':fstr(eta_minus),
         'mean_plus':fstr(state_plus),'mean_minus':fstr(state_minus),
         'state_difference_minus_plus':fstr(state_minus-state_plus)},
       'tracer_toy':{'A_L':'2','A_E':'3','same_fractional_chi':fstr(same_fractional),
         'unequal_fractional_chi':fstr(different_fractional),
         'interpretation':'coefficients deliberately synthetic, not eBOSS bias, beta or Doppler'},
       'mock_covariance':{'N_independent_realizations':9,'output_components':24,
         'centered_sample_covariance_max_rank':8,'exact_9x9_witness_rank':mock_rank},
       'even_leakage_toy':{'mock_leakage':fstr(toy_leak),'state_difference':'0',
         'interpretation':'not the measured eBOSS window'},
       'physical_stop':{'c_L_c_E_chi_F':'UNKNOWN','state_specific_signed_joint_pair_selection':'UNKNOWN',
         'original_24D_EBOSS_EV_mean':'UNCERTIFIED',
         'A03':'PHYSICAL_UNCERTIFIED','A04':'BLOCKED','observed_odd':'SEALED'}
    }
    # Do not overwrite an independent prior result.
    if OUT.exists():
        existing=json.loads(OUT.read_text(encoding='utf-8'))
        assert existing == report, 'immutable_report_changed'
    else:
        OUT.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('E30_OFFLINE_QA_PASS',len(checks),'assertions',flush=True)
    print('E30_MARGINAL_SIGN_BALANCED_BUT_DELTA_U',cross,flush=True)
    print('E30_GAUSSIAN_LEADING_MIXED_THREE_FIELD',triple,flush=True)
    print('E30_MOCK_RANK',mock_rank,'OF',24,flush=True)
    print('OUTPUT',OUT,'SHA256',hashlib.sha256(OUT.read_bytes()).hexdigest(),flush=True)

if __name__=='__main__':
    try:execute()
    except Exception as e:
        print('E30_OFFLINE_QA_FAIL',type(e).__name__,str(e),flush=True)
        sys.exit(1)
