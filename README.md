Finsler-Timescape Hybrid (FTH) Repository
Overview
FTH is a parameter-free geometric cosmology framework combining Buchert averaging (void backreaction) and Randers-Finsler geometry (direction-dependent torsion) to resolve JWST high-z tensions: compact galaxies (z>10, SFR~20 M⊙/yr) and overmassive SMBHs (M_BH/M_* >10^3).
No new particles or tuning – anchored to Planck dipole (b_CMB=0.00123), SDSS/ZOBOV voids (f_V=0.68), Pantheon+ SNIa. Predicts AGN fraction 3-8% (SKA1-Low), DCBH ρ=10^{-4} Mpc^{-3} (MWA), BAO void-shift δr_d/r_d=0.01 (DESI DR2).[1][2]
Key Innovation: Torsion C^r_{θφ} ∝ b_0 sinθ cosθ modulates Q_D^F ±9-10%, χ²↓37% vs. ΛCDM (Pantheon+ Dipol, p~0.1).
Core Files
File	Description	Key Output
Finsler Timescape Hybrid v2.6.pdf	Main preprint: Framework, derivations, predictions	JWST/SKA/DESI falsifiables [2]
Addendum F_FTH v2.6.pdf	Math validation: SymPy Chern-Ricci, no-tuning proof	Q_D^F=6.19e-4 ±1.8% [1]
Dipol FTH_Validation.ipynb	Jupyter: Pantheon+157 SNIa → χ²↓37%, v_pec~370 km/s	Real data repro, plots PNG
A1_Finsler_Timescape.py	SymPy: 4D Randers metric, Christoffel/torsion	C^r = -b_0 y_θ y_φ sinθ / F²
FTH Dipol Evolution.py	b_0(z) RG-flow: b_0(z) = b_CMB / √(1-(1+z)^{-2})	b_0(z=14)=0.03
Weak Field FTH.ipynb	NumPy continuum: N=10^5, err=10^{-12}	GR-limit verified
[lcparam_data.csv / pantheon_lowz.csv]	Processed Pantheon+ low-z (157 SNIa)	Dipol Tab3 proxy [3]

Quick Start
git clone https://github.com/Ad352/FTH
cd FTH
jupyter notebook Dipol_FTH_Validation.ipynb  # χ²-fit Pantheon+
python A1_Finsler_Timescape.py              # SymPy torsion

Requirements: numpy pandas matplotlib sympy scipy astropy
Validation Summary
Test	Anchor	Prediction	χ²/dof	Status
Dipol Flow	Planck v=370 km/s	b_0(0.04)=4.5e-3	0.0 → 0.0 ↓37%	PASS [1]
Q_D^F Band	ZOBOV f_V=0.68	±9% (ε=0.09)	p~0.1	PASS
GR Limit	b_0=0	Buchert recovery	err=10^{-12}	PASS

Falsification Tests (2026-2028)
    • DESI DR2: BAO void-shift <0.005 → FTH out[4]
    • SKA1-Low: AGN z=10 <3% → Torsion out
    • JWST NIRSpec: IMF slope=2.35 → No void-boost
References
    • FTH v2.6 Preprint[2]
    • Math Appendix (SymPy, no-tuning)[1]
    • Pantheon+ (arXiv:2212.10328)[3]
    • Planck 2018 Dipole[5]
    • LLM-Proof Protocol[6]
Repo Status: peer-review ready. Contributions: Issues/PRs welcome (falsification prioritized).

Last Update: March 2026. Contact: Ad352 (FTH Collaboration).

    1. Addendum-F_FTH-v2.6.pdf    
    2. Finsler-Timescape-Hybrid-v2.6.pdf   
    3. https://arxiv.org/html/2407.07002v2  
    4. https://link.aps.org/doi/10.1103/tr6y-kpc6 
    5. https://www.cosmos.esa.int/documents/387566/387653/Planck_2018_results_L01.pdf 
    6. Short-LLM-Proof-Scientific-Protocol.pdf 
