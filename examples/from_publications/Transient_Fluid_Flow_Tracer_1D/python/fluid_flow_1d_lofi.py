# This is the retrieved script seen by MATLAB
from retriever import *

# Set port to call
PORT = 5000; 


# Access Functions
state_solve = lambda *x, **kwargs: call_remote_function('state_solve', PORT, *x, **kwargs)
J = lambda *x: call_remote_function('J', PORT, *x)
Jz = lambda *x: call_remote_function('Jz', PORT, *x)
J_uu_apply = lambda *x: call_remote_function('J_uu_apply', PORT, *x)
J_zz_apply = lambda *x: call_remote_function('J_zz_apply', PORT, *x)
apply_solution_operator_z_jacobian_transpose = lambda *x: call_remote_function('apply_solution_operator_z_jacobian_transpose', PORT, *x)
apply_solution_operator_z_jacobian = lambda *x: call_remote_function('apply_solution_operator_z_jacobian', PORT, *x)
apply_rs_hessian = lambda *x: call_remote_function('apply_rs_hessian', PORT, *x)
misfit_gradient = lambda *x: call_remote_function('misfit_gradient', PORT, *x)
apply_misfit_hessian = lambda *x: call_remote_function('apply_misfit_hessian', PORT, *x)

# Access Variables
num_steps = get_remote_variable('num_steps', PORT)
mesh_coordinates = get_remote_variable('mesh_coordinates', PORT)
K_mat = get_remote_variable('K_mat', PORT)
K_mat_one = get_remote_variable('K_mat_one', PORT)
M_one = get_remote_variable('M_one', PORT)
M = get_remote_variable('M', PORT)

