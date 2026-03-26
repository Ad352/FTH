import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# SymPy for b0 evolution
z_sym, bCMB = sp.symbols('z bCMB')
b0_sym = bCMB * sp.sqrt(1 - sp.exp(-2*sp.ln(1+z_sym)))
print("SymPy b0(z=14):", sp.N(b0_sym.subs({z_sym:14, bCMB:0.00123}), 4))

# Numerical plot data
z = np.logspace(-1, 2, 100)
a = 1 / (1 + z)
bCMB_num = 0.00123
b0_num = bCMB_num * np.sqrt(1 - np.exp(-2 * np.log(1 + z)))

# Plot
fig, ax = plt.subplots(figsize=(8,5))
ax.semilogx(z, b0_num * 1000, 'b-', linewidth=2, label='b₀(z) × 10³')
ax.axhline(bCMB_num * 1000, color='r', linestyle='--', label='b_CMB (IR-Fixed-Point)')
ax.set_xlabel('Redshift z')
ax.set_ylabel('Dipol-Amplitude b₀ × 10³')
ax.set_title('FTH: Evolution des makroskopischen Dipols (Ricci-Fluss)')
ax.grid(True, alpha=0.3)
ax.legend()
plt.savefig('output/FTH_dipol_evolution.png', dpi=150, bbox_inches='tight')
plt.close()
print("Plot saved as FTH_dipol_evolution.png")