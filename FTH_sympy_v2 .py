import sympy as sp

# Coordinates & Symbols
t, r, th, ph, b0, eps, lam, bCMB, ell = sp.symbols('t r theta phi b0 epsilon lambda bCMB ell', real=True, positive=True)
yt, yr, yth, yph = sp.symbols('yt yr ytheta yph', real=True)
M_sym = sp.symbols('M', positive=True)
H0, MPl, T_Riem, a, QD, Omegam = sp.symbols('H0 M_Pl T_Riem a QD Omega_m', positive=True)

# 1. Randers-Finsler (spherical approx, +--- signature)
alpha2 = yt**2 - yr**2 - r**2 * yth**2 - r**2 * sp.sin(th)**2 * yph**2
alpha = sp.sqrt(alpha2)
beta = b0 * sp.cos(th) * yr
F = alpha + beta
F2 = F**2

# Finsler Metric components
g_rr = sp.diff(F2, yr, yr) / 2
g_tt = sp.diff(F2, yt, yt) / 2

# Christoffel schematic ^r_rr
Gamma_rr_r = sp.diff(g_rr, r) / (2 * g_rr)

# 2. Riemann Reg. C.6
R_rr = b0 / (3 * r**2)
R_reg = R_rr * sp.exp(-r**2 / (2 * eps**2))

# 3. RG-Flow S5.3
b0_ell = sp.Function('b0')(ell)
RG_eq = sp.Eq(b0_ell.diff(ell), lam * b0_ell * (b0_ell - bCMB))
RG_sol = sp.dsolve(RG_eq, b0_ell)

# 4. λ Ricci Proxy: Tr(R^i_jkb b^k) ~ H0 M_Pl / T_Riem * (b0 / r^2)
lambda_ricci = (H0 * MPl / T_Riem) * (b0 / r**2)

# 5. Torsion Torque C^r_θφ
C_r_thph = -b0 * sp.sin(th)
torque = C_r_thph * M_sym * r * (1 + b0)

# 6. Proper Time 2.3
H_D = sp.sqrt(Omegam / a**3 + QD / (3 * a**2))
dt_da = 1 / (a * H_D * (1 + QD/3))
t_void = sp.Integral(dt_da, (a, sp.Rational(1,15), 1))

# Outputs (LaTeX)
print('RG Equation:', sp.latex(RG_eq))
print('RG Solution:', sp.latex(RG_sol))
print('R_reg:', sp.latex(R_reg))
print('lambda Ricci:', sp.latex(lambda_ricci))
print('Torque:', sp.latex(torque))
print('t_void:', sp.latex(t_void))
print('g_rr sample:', sp.latex(g_rr.simplify()))

# Numerical subs example
subs = {b0:0.1, eps:1, r:0.001, lam:0.01, bCMB:0.0012, H0:70, MPl:1.22e19, T_Riem:1e-3, M_sym:1e5}
print('\\nNum lambda_ricci:', float(lambda_ricci.subs(subs)))
print('Num R_reg(r=0.001):', float(R_reg.subs({b0:0.1, eps:1, r:0.001})))
