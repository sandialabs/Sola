# IF PKG_CONFIG NOT FOUND
# try:
#     import matlab_env_fix
#     import os
#     os.environ = matlab_env_fix.data
#     import sys
#     sys.setdlopenflags(10)
# except ModuleNotFoundError:
#     pass

from fenics_helpers import * 
from pathlib import Path
from pyadjoint.reduced_functional_numpy import ReducedFunctionalNumPy
from scipy.io import loadmat
from scipy.linalg import block_diag


root_path = Path(__file__).parent
# Note that fenics_helpers assumes interval of [0, 1] to remove bug from ds

# Mesh Setup
N = 30;
unit_mesh = UnitIntervalMesh(N)
ds = generate_ds(unit_mesh)

# Set Function Spaces
P1 = FiniteElement('CG', unit_mesh.ufl_cell(), 1)
P2 = FiniteElement('CG', unit_mesh.ufl_cell(), 2)
U = FunctionSpace(unit_mesh, P2) # Velocity
K = FunctionSpace(unit_mesh, P1) # Tracer
mesh_coordinates = K.tabulate_dof_coordinates()

# Define Temporal Mesh
T = 0.1
num_steps = 25
total_dof = num_steps * (N+1)
dt = Constant(T/num_steps)
t = Constant(0);
gamma = Constant(0.05)
# reac_fn = lambda c: Constant(2) * (c+Constant(1))**2
reac_fn = lambda c: Constant(1) * c

# Retreive velocity from Timeseries
u_timeseries = TimeSeries(f"{root_path}/../../data/velocity_timeseries_midfi_1d")

# Store Mass Matrix for Future
M_one = assemble(TrialFunction(K) * TestFunction(K) * dx).array()
M = block_diag(*[M_one] * num_steps)
K_mat_one = assemble(TrialFunction(K).dx(0) * TestFunction(K).dx(0) * dx).array()
K_mat = block_diag(*[K_mat_one] * num_steps)


# Initial setup for inverse problem
k_terminal = fenics_convert(loadmat(f'{root_path}/../../data/terminal_state.mat', squeeze_me=True)["k_terminal"], "function", fun_space=K)
beta = Constant(1e-5)

# def state_solve(k0_input, return_type: Literal["vertex", "vector", "petsc", "function"], plot_k=False, annotate=True, verbose=False, return_all = False):
#     # Handle annotation and verbosity
#     if annotate: get_working_tape().clear_tape()
#     else: pause_annotation()
#     log_level = get_log_level()
#     set_log_active(verbose)

#     # Set Trial and Test Functions
#     k = TrialFunction(K)
#     v_k = TestFunction(K)

#     # Set initial conditions
#     k_list = [];
#     k0 = fenics_convert(k0_input, "function", K)
#     k_n = Function(K)
#     k_n.assign(k0)
#     if plot_k: plot(k_n)

#     # Solve the PDE with time-stepping
#     t.assign(0.0)

#     u_n = Function(U)
#     F_k = (1/dt * (k - k_n) * v_k + gamma*k.dx(0)*v_k.dx(0) + u_n*k.dx(0)*v_k + u_n.dx(0)*k*v_k + reac_fn(k)*v_k) * dx
#     a, L = lhs(F_k), rhs(F_k)

#     k_sol = Function(K)
#     for n in range(num_steps):
#         t.assign(float(t)+float(dt))
#         u_timeseries.retrieve(u_n.vector(), float(t))
#         solve(a == L, k_sol)
#         k_n.assign(k_sol)
#         if return_all: k_list.append(k_n.vector()[:])
#         if plot_k: plot(k_n)

#     # Revert annotation and verbosity
#     if not annotate: continue_annotation()
#     set_log_level(log_level)

#     # Return final state
#     if return_all: return np.array(k_list).flatten()
#     return fenics_convert(k_n, return_type, K)


def state_solve(k0_input, return_type: Literal["vertex", "vector", "petsc", "function"], plot_k=False, annotate=True, verbose=False, return_all = False):
    # Handle annotation and verbosity
    if annotate: get_working_tape().clear_tape()
    else: pause_annotation()
    log_level = get_log_level()
    set_log_active(verbose)

    # Set Trial and Test Functions
    k = Function(K)
    v_k = TestFunction(K)

    # Set initial conditions
    k_list = [];
    k0 = fenics_convert(k0_input, "function", K)
    k_n = Function(K)
    k_n.assign(k0)
    if plot_k: plot(k_n)

    # Solve the PDE with time-stepping
    t.assign(0.0)
    for n in range(num_steps):
        t.assign(float(t)+float(dt))
        u_n = Function(U)
        u_timeseries.retrieve(u_n.vector(), float(t))
        F_k = (1/dt * (k - k_n) * v_k + gamma*k.dx(0)*v_k.dx(0) + u_n*k.dx(0)*v_k + u_n.dx(0)*k*v_k + reac_fn(k)*v_k) * dx
        solve(F_k == 0, k, J=derivative(F_k, k))
        k_n.assign(k)
        if return_all: k_list.append(k_n.vector()[:])
        if plot_k: plot(k_n)

    # Revert annotation and verbosity
    if not annotate: continue_annotation()
    set_log_level(log_level)

    # Return final state
    if return_all: return np.array(k_list).flatten()
    return fenics_convert(k_n, return_type, K)

def state_solve_all_obj(k0, kt_in, annotate=True, verbose=False):
    # Handle annotation and verbosity
    if annotate: get_working_tape().clear_tape()
    else: pause_annotation()
    log_level = get_log_level()
    set_log_active(verbose)

    # Set Trial and Test Functions
    k = Function(K)
    v_k = TestFunction(K)

    # Set initial conditions
    J_output = 0.0;
    kt_in = np.linalg.solve(M_one, fenics_convert(kt_in, "vector").reshape(num_steps, N+1).T)
    k0 = fenics_convert(k0, "function", K)
    k_n = Function(K)
    k_n.assign(k0)

    # Solve the PDE with time-stepping
    t.assign(0.0)
    for n in range(num_steps):
        t.assign(float(t)+float(dt))
        u_n = Function(U)
        u_timeseries.retrieve(u_n.vector(), float(t))
        F_k = (1/dt * (k - k_n) * v_k + gamma*k.dx(0)*v_k.dx(0) + u_n*k.dx(0)*v_k + u_n.dx(0)*k*v_k + reac_fn(k)*v_k) * dx
        solve(F_k == 0, k, J=derivative(F_k, k))
        k_n.assign(k)
        J_output += assemble(inner(k_n, fenics_convert(kt_in[:, n], "function", K)) * dx)

    # Revert annotation and verbosity
    if not annotate: continue_annotation()
    set_log_level(log_level)

    # Return sum of inner product objectives with test vectors
    return J_output

def J(k0, kt):
    # Convert inputs to functions
    k0 = fenics_convert(k0, "function", fun_space=K)
    kt = fenics_convert(kt, "function", fun_space=K)
    val = 1.e4*(assemble(0.5*inner(kt - k_terminal, kt - k_terminal)*dx + 0.5 * beta * inner(k0.dx(0), k0.dx(0)) * dx))
    return val

def Jz(k0, kt):
    # Convert inputs to vectors
    k0 = fenics_convert(k0, "vector", fun_space=K)
    val = 1.e4*(float(beta) * K_mat_one @ k0)
    return val

J_hat_np = None;
def reduced_functional_J_hat(k0):
    global J_hat_np;
    # Set global beta & k_terminal
    k0 = fenics_convert(k0, "function", fun_space=K)
    kt = state_solve(k0, "function")

    J_inv = J(k0, kt)
    control = Control(k0)
    J_hat = ReducedFunctional(J_inv, control)
    J_hat_np = ReducedFunctionalNumPy(J_hat)
    return J_hat_np

def callback_call(J_hat, full=True):
    c_test = Function(K)

    def callback_fn(x):
        if callback_fn.iteration == 0:
            # Print the header row
            print(f"{'Iteration':<10} {'Objective Value':<20} {'Gradient Norm':<15}")
            print("-" * 45)
        callback_fn.iteration += 1
        if full:
            c_test.vector()[:] = x.flatten()
            print(f"{callback_fn.iteration:<10} {J_hat(c_test):<20.6f} {norm(J_hat.derivative()):<15.6f}")
        else:
            print(f"{callback_fn.iteration:<10}")
    callback_fn.iteration = 0
    
    return callback_fn

# Functions to be called by Matlab Interface
# ------------------------------------------

def apply_solution_operator_z_jacobian_transpose(kt_in, k0):
    # Multiply by inv(M) to compensate for assembly in L2
    k0 = fenics_convert(k0, "function", fun_space=K)
    return compute_gradient(state_solve_all_obj(k0, kt_in), Control(k0)).vector()[:]

def misfit_gradient(kt, k0):
    # Backwards Compatible for more
    kt = fenics_convert(kt, "vector", K)
    if kt.size == N + 1: return 1.e4 * M_one @ (kt - fenics_convert(k_terminal, "vector"));
    unpadded_vec = 1.e4 * M_one @ (kt[-N-1:] - fenics_convert(k_terminal, "vector"));
    return np.pad(unpadded_vec, (kt.size-N-1, 0), 'constant');

def apply_misfit_hessian(kt_in, k0, kt):
    # Backwards Compatible for more
    kt_in = fenics_convert(kt_in, "vector", K)
    if kt_in.size == N + 1: raise Exception("Incorrect kt_in size!")
    unpadded_vec = 1.e4 * M_one @ fenics_convert(kt_in[-N-1:], "vector", K);
    return np.pad(unpadded_vec, (kt_in.size-N-1, 0), 'constant');

def apply_rs_hessian(k0_in_mat, k0):
    def apply_rs_hessian_single(k0_in):
        k0_in = fenics_convert(k0_in, "function", fun_space=K)
        kt_in = state_solve(k0_in, return_type="function", return_all=False)
        kt_copy = Function(K)
        kt_copy.vector()[:] = kt_in.vector()[:]
        tmp = assemble(inner(kt_in, kt_copy) * dx)
        tmp2 = 1e4*compute_gradient(tmp, Control(k0_in)).vector()[:] + Jz(k0_in, None)
        return tmp2

    if k0_in_mat.ndim == 1 or k0_in_mat.shape[1] == 1: return apply_rs_hessian_single(k0_in_mat)
    return np.apply_along_axis(apply_rs_hessian_single, 0, k0_in_mat)


def apply_solution_operator_z_jacobian(k0_in, k0):
    return state_solve(k0_in, return_type="vector", return_all=True)
   

print(f"Succesfully imported {__name__}")

