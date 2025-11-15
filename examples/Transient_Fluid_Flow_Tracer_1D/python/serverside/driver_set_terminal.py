from pathlib import Path
from scipy.io import savemat
import os

root_path = Path(__file__).parent

if os.path.exists(f"{root_path}/../../data/velocity_timeseries_midfi_1d.h5"):
    os.remove(f"{root_path}/../../data/velocity_timeseries_midfi_1d.h5")

import fluid_flow_1d_hifi_eval as hifi

# Solve initial iterate for inverse problem (used for taping adjoints)
# k0_hifi = hifi.interpolate(hifi.Expression("x[0]*(1-x[0])*(9-10*x[0])", degree=1), hifi.C)
# k0_hifi = hifi.interpolate(hifi.Expression("10*pow(x[0],2)*pow(x[0]-1,2)", degree=1), hifi.C)
k0_hifi = hifi.interpolate(hifi.Expression("exp(-10*pow(x[0]-0.5,2))", degree=1), hifi.C)
k_terminal = hifi.state_solve(k0_hifi, return_type = "vector")
hifi.state_solve(0*k0_hifi.vector()[:], return_type = "vector", store_midfi=True); # 

# Save the array to a .mat file
savemat(f'{root_path}/../../data/terminal_state.mat', {'k_terminal': k_terminal}, oned_as='column')
savemat(f'{root_path}/../../data/hifi_optim_sol.mat', {'k0_hifi': hifi.fenics_convert(k0_hifi, "vector")}, oned_as='column')
print("Saved.")