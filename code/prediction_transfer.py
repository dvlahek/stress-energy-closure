#!/usr/bin/env python3
"""Model-level metric-response prediction for the stress-energy-matched pair.

Uses the same dimensionless validation convention as injectivity_tomography.py.
The plotted quantities are not observational sensitivity forecasts.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from injectivity_tomography import PMAX, TAUMAX, N_P, N_TAU, construct_matched_pair, metric_transfer


def main():
    out = Path("prediction_response_output")
    out.mkdir(exist_ok=True)
    p = np.linspace(0.0, PMAX, N_P)
    taus = np.linspace(0.0, TAUMAX, N_TAU)
    pair = construct_matched_pair(p, taus)
    omega, gp = metric_transfer(taus, pair["K_plus"])
    _, gm = metric_transfer(taus, pair["K_minus"])

    complex_sep = 100.0 * np.abs(gp-gm) / np.abs(gp)
    amp_sep = 100.0 * np.abs(np.abs(gp)-np.abs(gm)) / np.abs(gp)
    phase_deg = np.degrees(np.angle(gp/gm))

    np.savetxt(
        out / "prediction_transfer_response.csv",
        np.column_stack([omega, complex_sep, amp_sep, phase_deg]),
        delimiter=",",
        header="omega,relative_complex_response_difference_percent,relative_amplitude_difference_percent,phase_difference_degrees",
        comments="",
    )

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.figure(figsize=(7.0, 4.4))
    plt.plot(omega, complex_sep, label="complex response")
    plt.plot(omega, amp_sep, label="amplitude")
    plt.xlabel(r"$\omega$")
    plt.ylabel("relative response difference (%)")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out / "prediction_transfer_response.pdf")
    plt.savefig(out / "prediction_transfer_response.png", dpi=300)
    plt.close()

    print(f"max complex response difference: {complex_sep.max():.6f}%")
    print(f"max amplitude difference: {amp_sep.max():.6f}%")
    print(f"max |phase difference|: {np.max(np.abs(phase_deg)):.6f} deg")


if __name__ == "__main__":
    main()
