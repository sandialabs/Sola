import fluid_flow_1d_lofi as lofi
from scipy.io import savemat

# Solve initial iterate for inverse problem (used for taping adjoints)
k0_guess = lofi.interpolate(lofi.Expression("2*(0.4 < x[0] && x[0] < 0.6)", degree=1), lofi.K)
k_n = lofi.state_solve(k0_guess, return_type = "function", plot_k=False, verbose=False, annotate=True);

# Set up inverse problem
lofi.beta.assign(1e-5)
J_inv = lofi.J(k0_guess, k_n)
control = lofi.Control(k0_guess)
J_hat = lofi.ReducedFunctional(J_inv, control)

# Solve Inverse Problem
with lofi.stop_verbose():
    k0_opt_lofi = lofi.minimize(J_hat, method="Newton-CG", callback=lofi.callback_call(J_hat), tol=1e-3, options={"disp": True})

k_opt_lofi = lofi.state_solve(k0_opt_lofi, return_type = "vector", plot_k=False, verbose=False, annotate=True);

# Save the array to a .mat file
savemat('data/lofi_optim_sol.mat', {'k0_opt_lofi': lofi.fenics_convert(k0_opt_lofi, "vector"), 'k_opt_lofi': k_opt_lofi}, oned_as='column')