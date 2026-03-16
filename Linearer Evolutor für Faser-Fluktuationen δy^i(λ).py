import sympy as sp
from scipy.integrate import solve_ivp
import numpy as np

print("=== Python Modul: Linearer Evolutor für Faser-Fluktuationen δy^i(λ) ===")

# SymPy: Linear Matrix A^i_j
lambda_sym = sp.symbols('lambda')
delta_y1, delta_y2, delta_y3 = sp.symbols('delta_y1 delta_y2 delta_y3')
C_r_th_ph_sym = sp.symbols('C_r_theta_phi')
v1, v2, v3 = sp.symbols('v1 v2 v3')

A_matrix = sp.Matrix([[-C_r_th_ph_sym * v3, 0, -C_r_th_ph_sym * v1],
                      [0, 0, 0],
                      [C_r_th_ph_sym * v2, -C_r_th_ph_sym * v1, 0]])
print("A-Matrix:", sp.latex(A_matrix))

# NumPy ODE
def fiber_fluct_ode(t, delta_y, C_params, v):
    """dδy/dλ = A δy mit A^i_j = -C^i_jk v^k"""
    b0, sin_th, F_inv2 = C_params
    C_rthph = -b0 * sin_th * F_inv2
    A = np.array([[-C_rthph * v[2], 0, -C_rthph * v[0]],
                  [0, 0, 0],
                  [C_rthph * v[1], -C_rthph * v[0], 0]])
    return A @ delta_y

# Demo-Integration
C_params = (0.03, np.sin(1.0), 1.0)  # b0, sinθ, 1/F^2
v_demo = np.array([1.0, 0.5, 0.2])
y0 = [1e-3, 5e-4, 2e-4]

sol = solve_ivp(fiber_fluct_ode, [0, 10], y0, args=(C_params, v_demo),
                method='RK45', rtol=1e-8, dense_output=True)

print("δy(λ=10):", sol.y[:,-1])
print("Norm-Wachstum:", np.linalg.norm(sol.y[:,-1]) / np.linalg.norm(y0))
reynolds = np.trace(np.outer(sol.y[:,-1], sol.y[:,-1]))
print("Reynolds ~ ν_eff:", reynolds)

print("Modul ready für FTH Appendix A.2: PDE → ODE → Reynolds-Stress.")