Finsler-Timescape Hybrid (FTH) Repository
Overview
FTH is a parameter-free geometric cosmology framework combining Buchert averaging (void backreaction) and Randers-Finsler geometry (direction-dependent torsion) to resolve JWST high-z tensions: compact galaxies (z>10, SFR~20 M⊙/yr) and overmassive SMBHs (M_BH/M_* >10^3).

FTH is a parameter-free geometric cosmology framework combining Buchert averaging (void backreaction) and Randers-Finsler geometry (direction-dependent torsion) to resolve JWST high-z tensions: compact galaxies (z>10, SFR>20 M⊙/yr) and overmassive SMBHs
Finsler-Timescape Hybrid (FTH) 

Repository Overview

Core Files
File	Description	Key Output

Finsler Timescape Hybrid v2.6.pdf	Main preprint: Framework, derivations, predictions	JWST/SKA/DESI falsifiables

Addendum F_FTH v2.6.pdf	Math validation: SymPy Chern-Ricci, no-tuning proof	Q_D^F=6.19e-4 ±1.8%

Dipol FTH_Validation.ipynb	Jupyter: Pantheon+157 SNIa → χ²↓37%, v_pec~370 km/s	Real data repro, plots PNG

A1_Finsler_Timescape.py	SymPy: 4D Randers metric, Christoffel/torsion	C^r = -b_0 y_θ y_φ sinθ / F²

FTH Dipol Evolution.py	b_0(z) RG-flow: b_0(z) = b_CMB / √(1-(1+z)^{-2})	b_0(z=14)=0.03

Weak Field FTH.ipynb	NumPy continuum: N=10^5, err=10^{-12}	GR-limit verified

[lcparam_data.csv / pantheon_lowz.csv]	Processed Pantheon+ low-z (157 SNIa)	Dipol Tab3 proxy


