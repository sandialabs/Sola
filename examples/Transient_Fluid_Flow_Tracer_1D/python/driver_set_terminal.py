import fluid_flow_1d_hifi_eval as hifi
from pathlib import Path
from scipy.io import savemat

root_path = Path(__file__).parent


# Solve initial iterate for inverse problem (used for taping adjoints)
k0_hifi = hifi.interpolate(hifi.Expression("x[0]*(1-x[0])*(9-10*x[0])", degree=1), hifi.C)
k_terminal = hifi.state_solve(k0_hifi, return_type = "vector")

# Save the array to a .mat file
savemat(f'{root_path}/../data/terminal_state.mat', {'k_terminal': k_terminal}, oned_as='column')
savemat(f'{root_path}/../data/hifi_optim_sol.mat', {'k0_hifi': hifi.fenics_convert(k0_hifi, "vector")}, oned_as='column')