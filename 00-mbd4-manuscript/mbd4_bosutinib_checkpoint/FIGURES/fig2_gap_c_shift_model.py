"""
GAP C — IN SILICO: Checkpoint Depletion Threshold Shift Model
==============================================================
We CAN model this computationally. Here's the logic:

FROM THE MANUSCRIPT (real data, receipt-backed):
- MBD4-LOF shifts ATRi (ceralasertib) LN_IC50 by −0.74 log-units (p=0.034)
- After TP53 control: shift is −1.063 log-units (p=0.008)
- WEE1i (adavosertib) shift: −0.509 log-units (p=0.076)

WHAT THIS MEANS:
- A shift of −0.74 in LN_IC50 means the functional IC50 drops by a factor of
  exp(0.74) ≈ 2.1x
- A shift of −1.063 means the IC50 drops by exp(1.063) ≈ 2.9x

SO: If the enzymatic IC50 for bosutinib vs Wee1 is 644 nM, and MBD4-LOF
    shifts checkpoint inhibitor sensitivity by 2.1–2.9x, then the FUNCTIONAL
    threshold in MBD4-LOF cells would be:
    
    Conservative (0.5 log shift): 644 / exp(0.5)  = 644 / 1.65 = 390 nM
    Central     (0.74 log shift): 644 / exp(0.74) = 644 / 2.10 = 307 nM
    Aggressive  (1.06 log shift): 644 / exp(1.06) = 644 / 2.89 = 223 nM

THEN re-run the PK selectivity model with these shifted thresholds.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

def run_simulation(dose_mg, median_cmin_nM, t_tumor_mean, t_marrow_mean,
                   t_tumor_sd, t_marrow_sd, cv=0.30, n=10000):
    sigma_ln = np.sqrt(np.log(cv**2 + 1))
    mu_ln = np.log(median_cmin_nM)
    cmin = np.random.lognormal(mean=mu_ln, sigma=sigma_ln, size=n)
    
    t_tumor = np.random.normal(t_tumor_mean, t_tumor_sd, n)
    t_marrow = np.random.normal(t_marrow_mean, t_marrow_sd, n)
    t_marrow = np.maximum(t_marrow, 500)
    
    in_kill_zone = (cmin > t_tumor) & (cmin <= t_marrow)
    marrow_tox = cmin > t_marrow
    sub_therapeutic = cmin <= t_tumor
    
    return {
        'dose': dose_mg,
        'median_cmin': median_cmin_nM,
        'pct_sub': np.mean(sub_therapeutic) * 100,
        'pct_kill': np.mean(in_kill_zone) * 100,
        'pct_tox': np.mean(marrow_tox) * 100,
    }


def main():
    np.random.seed(42)
    N = 10000
    
    # From PMC4111673 (Beeharry et al.)
    BOSUTINIB_WEE1_IC50 = 644   # nM (enzymatic)
    BOSUTINIB_CHK1_IC50 = 785   # nM (enzymatic)
    
    # Use the lower of the two (Wee1) as the functional checkpoint threshold
    ENZYMATIC_IC50 = BOSUTINIB_WEE1_IC50
    
    # From manuscript GDSC2 data:
    # MBD4-LOF shifts checkpoint inhibitor LN_IC50 by these amounts
    shifts = {
        'WT (no shift)':           0.0,
        'WEE1i trend (Δ=−0.51)':  0.509,
        'ATRi signal (Δ=−0.74)':  0.738,
        'TP53-ctrl (Δ=−1.06)':    1.063,
    }
    
    # Marrow threshold stays at 1500 nM (CD34+ HSPC data, drug-agnostic)
    T_MARROW = 1500
    T_MARROW_SD = 300
    
    # Dose-Cmin mapping (validated Pfizer PK)
    doses = [300, 400, 500, 600, 700]
    cmin_medians = [167, 222, 277, 333, 389]
    
    print("=" * 80)
    print("GAP C — IN SILICO: MBD4-LOF Checkpoint Depletion Shift Model")
    print("=" * 80)
    print(f"Enzymatic IC50 (Bosutinib→Wee1): {ENZYMATIC_IC50} nM (PMC4111673)")
    print(f"Marrow Threshold: {T_MARROW} nM (CD34+ HSPC IC50 proxy)")
    print()
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    for idx, (scenario_name, shift_magnitude) in enumerate(shifts.items()):
        # Apply the shift: functional_threshold = enzymatic / exp(shift)
        if shift_magnitude > 0:
            functional_threshold = ENZYMATIC_IC50 / np.exp(shift_magnitude)
        else:
            functional_threshold = ENZYMATIC_IC50
        
        functional_sd = functional_threshold * 0.20  # 20% CV on threshold
        
        print(f"--- Scenario: {scenario_name} ---")
        print(f"  Shift: −{shift_magnitude:.3f} log-units → exp({shift_magnitude:.3f}) = {np.exp(shift_magnitude):.2f}x")
        print(f"  Functional T_tumor: {ENZYMATIC_IC50} / {np.exp(shift_magnitude):.2f} = {functional_threshold:.0f} nM")
        print(f"  {'Dose':>6} | {'Cmin':>8} | {'Sub-Ther':>10} | {'Kill Zone':>10} | {'Marrow Tox':>10}")
        print(f"  {'-'*60}")
        
        results = []
        for dose, cmin_med in zip(doses, cmin_medians):
            res = run_simulation(
                dose, cmin_med,
                functional_threshold, T_MARROW,
                functional_sd, T_MARROW_SD,
                n=N
            )
            results.append(res)
            flag = " ✓" if res['pct_kill'] > 50 else ""
            print(f"  {dose:>5}mg | {cmin_med:>6} nM | {res['pct_sub']:>9.1f}% | {res['pct_kill']:>9.1f}% | {res['pct_tox']:>9.1f}%{flag}")
        print()
        
        # Plot
        ax = axes[idx]
        d = [r['dose'] for r in results]
        k = [r['pct_kill'] for r in results]
        t = [r['pct_tox'] for r in results]
        
        ax.plot(d, k, 'go-', linewidth=2, markersize=8, label='Kill Zone %')
        ax.plot(d, t, 'r^-', linewidth=2, markersize=8, label='Marrow Tox %')
        ax.axhline(50, color='green', linestyle=':', alpha=0.5, label='50% efficacy')
        ax.axhline(5, color='red', linestyle=':', alpha=0.5, label='5% tox limit')
        ax.set_title(f'{scenario_name}\nT_tumor = {functional_threshold:.0f} nM', fontsize=11, fontweight='bold')
        ax.set_xlabel('Dose (mg)')
        ax.set_ylabel('% Patients')
        ax.set_ylim(-2, 100)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.suptitle('GAP C In Silico: How MBD4-LOF Checkpoint Depletion Shifts the Kill Zone\n'
                 '(Shift factors from GDSC2 ceralasertib/adavosertib data applied to bosutinib Wee1 IC50)',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    
    out_dir = '/Users/fahadkiani/.gemini/antigravity/brain/a1668b4b-113c-4e7b-a0cb-74316bcfff82'
    plot_path = os.path.join(out_dir, 'gap_c_in_silico.png')
    plt.savefig(plot_path, dpi=150)
    print(f"Plot saved: {plot_path}")
    
    # Summary
    print("=" * 80)
    print("BOTTOM LINE")
    print("=" * 80)
    print()
    print("At standard 500mg (Cmin=277nM):")
    for scenario_name, shift_magnitude in shifts.items():
        ft = ENZYMATIC_IC50 / np.exp(shift_magnitude) if shift_magnitude > 0 else ENZYMATIC_IC50
        res = run_simulation(500, 277, ft, T_MARROW, ft*0.20, T_MARROW_SD, n=N)
        verdict = "VIABLE" if res['pct_kill'] > 30 else "MARGINAL" if res['pct_kill'] > 10 else "NOT VIABLE"
        print(f"  {scenario_name:30s}: T_tumor={ft:>5.0f}nM → Kill={res['pct_kill']:.1f}%, Tox={res['pct_tox']:.1f}% [{verdict}]")
    
    print()
    print("INTERPRETATION:")
    print("  - Without MBD4-LOF shift (WT): bosutinib cannot reach Wee1 IC50. Dead hypothesis.")
    print("  - With WEE1i-class shift (0.51 log): T_tumor drops to ~390 nM. Still marginal at 500mg.")
    print("  - With ATRi-class shift (0.74 log): T_tumor drops to ~307 nM. Borderline at 500mg.")  
    print("  - With TP53-controlled shift (1.06 log): T_tumor drops to ~223 nM. Viable at 500mg.")
    print()
    print("  The hypothesis LIVES if and only if MBD4-LOF creates a >0.7 log-unit shift")
    print("  for bosutinib's checkpoint activity — same magnitude as ceralasertib in GDSC2.")
    print("  This is plausible but unproven. Wet lab required to confirm.")


if __name__ == '__main__':
    main()
