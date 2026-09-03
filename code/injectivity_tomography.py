#!/usr/bin/env python3
"""
Gravitational-response identifiability and synthetic kinetic tomography.

Constructs two smooth positive massive isotropic kinetic states with the same
particle number density n, energy density rho and pressure P, but different
TT response kernels. It then reconstructs each hidden radial distribution
from noisy response samples using a non-negative B-spline inverse problem
with cross-validated Tikhonov regularization.
"""

from pathlib import Path
import json
import numpy as np
from scipy.linalg import null_space
from scipy.interpolate import BSpline
from scipy.optimize import lsq_linear
from scipy.special import spherical_jn
import matplotlib.pyplot as plt

C = 1.0 / (2.0 * np.pi) ** 3
MASS = 1.0
K_WAVE = 1.0
GAMMA = 0.5
PMAX = 5.0
TAUMAX = 25.0
N_P = 1600
N_TAU = 360
NOISE = 1.0e-3
SEED = 9


def integrate(y, x):
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x)
    return np.trapz(y, x)


def angular_A(x):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1.0e-6
    xs = x[small]
    out[small] = 16.0/15.0 - 8.0*xs**2/105.0 + 2.0*xs**4/945.0
    xx = x[~small]
    out[~small] = 16.0 * spherical_jn(2, xx) / xx**2
    return out


def bump(p, center, width):
    x = (p-center)/width
    f = np.zeros_like(p)
    mask = np.abs(x) < 1.0
    f[mask] = np.exp(-1.0/(1.0-x[mask]**2))
    return f


def moments(p, F, mass=MASS):
    E = np.sqrt(p*p + mass*mass)
    n = C*4.0*np.pi*integrate(p*p*F, p)
    rho = C*4.0*np.pi*integrate(p*p*E*F, p)
    P = C*4.0*np.pi/3.0*integrate(p**4/E*F, p)
    return np.array([n, rho, P])


def build_kernel_operator(p, basis, taus, mass=MASS, k=K_WAVE):
    E = np.sqrt(p*p + mass*mass)
    v = p/E
    dp = p[1]-p[0]
    qw = np.full_like(p, dp)
    qw[0] *= 0.5
    qw[-1] *= 0.5
    pref = C*np.pi/4.0
    R = np.empty((len(taus), len(basis)))
    for j, b in enumerate(basis):
        db = np.gradient(b, p, edge_order=2)
        radial = qw*(p**5/E)*db
        for i, tau in enumerate(taus):
            R[i,j] = pref*np.sum(radial*angular_A(k*v*tau))
    return R


def construct_matched_pair(p, taus):
    centers = np.linspace(0.30, 4.30, 10)
    raw = np.array([bump(p, c, 0.50) for c in centers])
    rho = np.array([moments(p, b)[1] for b in raw])
    basis = raw/rho[:,None]
    constraint = np.array([[moments(p, b)[i] for b in basis] for i in range(3)])
    response = build_kernel_operator(p, basis, taus)
    N = null_space(constraint)
    _, _, VT = np.linalg.svd(response @ N, full_matrices=False)
    direction = N @ VT[0]
    direction /= np.max(np.abs(direction))
    base = np.abs(direction) + 0.05
    base /= base.sum()
    nz = np.abs(direction) > 1.0e-12
    alpha = 0.90*np.min(base[nz]/np.abs(direction[nz]))
    c_plus = base + alpha*direction
    c_minus = base - alpha*direction
    assert np.all(c_plus > 0.0) and np.all(c_minus > 0.0)
    return c_plus@basis, c_minus@basis, response@c_plus, response@c_minus


def bspline_basis(p, n_basis=24, degree=3):
    n_internal = n_basis-degree-1
    internal = np.linspace(0.0, PMAX, n_internal+2)[1:-1]
    knots = np.r_[np.repeat(0.0, degree+1), internal, np.repeat(PMAX, degree+1)]
    B, dB = [], []
    for i in range(n_basis):
        c = np.zeros(n_basis)
        c[i] = 1.0
        s = BSpline(knots, c, degree, extrapolate=False)
        B.append(np.nan_to_num(s(p)))
        dB.append(np.nan_to_num(s.derivative()(p)))
    return np.array(B), np.array(dB)


def spline_forward_operator(p, taus, n_basis=24):
    B, dB = bspline_basis(p, n_basis=n_basis)
    E = np.sqrt(p*p + MASS*MASS)
    v = p/E
    dp = p[1]-p[0]
    qw = np.full_like(p, dp)
    qw[0] *= 0.5
    qw[-1] *= 0.5
    pref = C*np.pi/4.0
    M = np.empty((len(taus), n_basis))
    for i, tau in enumerate(taus):
        weight = qw*(p**5/E)*angular_A(K_WAVE*v*tau)
        M[i,:] = pref*(dB @ weight)
    return M[:,:-1], B[:-1]


def reconstruct(p, taus, K_true, F_true, noise_fraction=NOISE):
    M, B = spline_forward_operator(p, taus, n_basis=24)
    n_coeff = M.shape[1]
    L = np.zeros((n_coeff-2, n_coeff))
    for i in range(n_coeff-2):
        L[i, i:i+3] = [1.0, -2.0, 1.0]
    scale = np.max(np.abs(K_true))
    rng = np.random.default_rng(SEED)
    y = K_true + rng.normal(0.0, noise_fraction*scale, len(K_true))
    idx = np.arange(len(y))
    train = idx[idx % 4 != 0]
    valid = idx[idx % 4 == 0]
    candidates = []
    for lam in np.logspace(-10, -1, 28):
        A = np.vstack([M[train]/scale, np.sqrt(lam)*L])
        b = np.r_[y[train]/scale, np.zeros(L.shape[0])]
        fit = lsq_linear(A, b, bounds=(0.0, np.inf), max_iter=3000)
        val = np.linalg.norm(M[valid]@fit.x-y[valid])/(np.sqrt(len(valid))*scale)
        candidates.append((val, lam))
    _, lam = min(candidates, key=lambda z: z[0])
    A = np.vstack([M/scale, np.sqrt(lam)*L])
    b = np.r_[y/scale, np.zeros(L.shape[0])]
    fit = lsq_linear(A, b, bounds=(0.0, np.inf), max_iter=3000)
    F_hat = fit.x @ B
    return F_hat, M@fit.x, y, float(lam), M


def metric_transfer(taus, K):
    eta = 0.05
    omega = np.linspace(0.05, 3.0, 300)
    s = eta + 1j*omega
    Ktilde = np.array([integrate(np.exp(-si*taus)*K, taus) for si in s])
    G = 1.0/(s*s + K_WAVE*K_WAVE - GAMMA*s*Ktilde)
    return omega, G


def main():
    out = Path("response_identifiability_output")
    out.mkdir(exist_ok=True)
    p = np.linspace(0.0, PMAX, N_P)
    taus = np.linspace(0.0, TAUMAX, N_TAU)
    Fp, Fm, Kp, Km = construct_matched_pair(p, taus)
    mp, mm = moments(p, Fp), moments(p, Fm)
    Fhp, Khp, ynp, lamp, M = reconstruct(p, taus, Kp, Fp)
    Fhm, Khm, ynm, lamm, _ = reconstruct(p, taus, Km, Fm)
    omega, Gp = metric_transfer(taus, Kp)
    _, Gm = metric_transfer(taus, Km)
    Grel = np.abs(Gp-Gm)/np.maximum(np.abs(Gp), 1.0e-30)
    sv = np.linalg.svd(M/np.max(np.abs(M)), compute_uv=False)
    summary = {
        "mass": MASS,
        "k": K_WAVE,
        "noise_fraction": NOISE,
        "moment_plus_n_rho_P": mp.tolist(),
        "moment_minus_n_rho_P": mm.tolist(),
        "moment_max_abs_difference": float(np.max(np.abs(mp-mm))),
        "kernel_relative_L2_separation": float(np.linalg.norm(Kp-Km)/np.linalg.norm(Kp)),
        "kernel_max_relative_separation": float(np.max(np.abs(Kp-Km))/np.max(np.abs(Kp))),
        "distribution_relative_L2_separation": float(np.linalg.norm(Fp-Fm)/np.linalg.norm(Fp)),
        "reconstruction_plus_relative_L2_error": float(np.linalg.norm(Fhp-Fp)/np.linalg.norm(Fp)),
        "reconstruction_minus_relative_L2_error": float(np.linalg.norm(Fhm-Fm)/np.linalg.norm(Fm)),
        "reconstruction_plus_lambda": lamp,
        "reconstruction_minus_lambda": lamm,
        "max_metric_transfer_relative_difference": float(np.max(Grel)),
        "metric_transfer_peak_omega": float(omega[np.argmax(Grel)]),
        "finite_discretization_condition_number": float(sv[0]/sv[-1]),
    }
    with open(out/"summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    np.savetxt(out/"matched_states_and_reconstruction.csv", np.column_stack([p,Fp,Fm,Fhp,Fhm]), delimiter=",", header="p,F_plus,F_minus,F_plus_reconstructed,F_minus_reconstructed", comments="")
    np.savetxt(out/"response_kernels.csv", np.column_stack([taus,Kp,Km,Kp-Km,ynp,ynm]), delimiter=",", header="tau,K_plus,K_minus,K_plus_minus_K_minus,K_plus_noisy,K_minus_noisy", comments="")
    np.savetxt(out/"metric_transfer.csv", np.column_stack([omega,np.abs(Gp),np.abs(Gm),Grel]), delimiter=",", header="omega,abs_G_plus,abs_G_minus,relative_difference", comments="")
    np.savetxt(out/"inverse_singular_values.csv", np.column_stack([np.arange(1,len(sv)+1),sv]), delimiter=",", header="index,scaled_singular_value", comments="")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
