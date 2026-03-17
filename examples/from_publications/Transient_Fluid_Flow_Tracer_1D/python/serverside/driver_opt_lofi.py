import fluid_flow_1d_lofi_linear as lofi
import os
from scipy.io import savemat, loadmat

# Solve initial iterate for inverse problem (used for taping adjoints)
try: 
    # k0_guess = lofi.fenics_convert(loadmat(f'{lofi.root_path}/../../data/lofi_optim_sol.mat')['k0_opt_lofi'], 'function', lofi.K)
    k0_guess = lofi.fenics_convert(loadmat(f'{lofi.root_path}/../../data/hifi_optim_sol.mat')['k0_hifi'], 'function', lofi.K)
except:
    k0_guess = lofi.interpolate(lofi.Expression("2*(0.4 < x[0] && x[0] < 0.6)", degree=1), lofi.K)
k_n = lofi.state_solve(k0_guess, return_type = "function", plot_k=False, verbose=False, annotate=True);

# Set up inverse problem
J_inv = lofi.J(k0_guess, k_n)
control = lofi.Control(k0_guess)
J_hat = lofi.ReducedFunctional(J_inv, control)

# Solve Inverse Problem
with lofi.stop_verbose():
    k0_opt_lofi = lofi.minimize(J_hat, method="Newton-CG", callback=lofi.callback_call(J_hat), tol=1, options={"disp": True, "cg_tol_modifier": 1e-5})
    # Note: The cg_tol_modifier is custom (multiplies to eta in scipy.optimize._optimize.py line 2114)
    

k_opt_lofi = lofi.state_solve(k0_opt_lofi, return_type = "vector", plot_k=False, verbose=False, annotate=True);

# Save the array to a .mat file
savemat(f'{lofi.root_path}/../../data/lofi_optim_sol.mat', {'k0_opt_lofi': lofi.fenics_convert(k0_opt_lofi, "vector"), 'k_opt_lofi': k_opt_lofi}, oned_as='column')
print("Saved.")
os._exit(0)