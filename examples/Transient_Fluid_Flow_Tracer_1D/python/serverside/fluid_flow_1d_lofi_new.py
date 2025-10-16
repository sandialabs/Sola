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
from typing import Callable, Literal, Tuple

# -------------------------------------------------------------------------
# Global setup
# -------------------------------------------------------------------------

root_path = Path(__file__).parent
# Note that fenics_helpers assumes interval of [0, 1] to remove bug from ds

# Mesh Setup
N = 30
unit_mesh = UnitIntervalMesh(N)
ds = generate_ds(unit_mesh)

# Set Function Spaces
P1 = FiniteElement('CG', unit_mesh.ufl_cell(), 1)
P2 = FiniteElement('CG', unit_mesh.ufl_cell(), 2)
U = FunctionSpace(unit_mesh, P2)   # Velocity
K = FunctionSpace(unit_mesh, P1)   # Tracer
mesh_coordinates = K.tabulate_dof_coordinates()

# Temporal discretisation
T = 0.1
num_steps = 25
dt = Constant(T / num_steps)
t = Constant(0.0)
gamma = Constant(0.025)
beta = Constant(1e-5)
reac_fn = lambda c: Constant(10) * c     # - Constant(5)

# External data
u_timeseries = TimeSeries(f"{root_path}/../../data/velocity_timeseries_midfi_1d")
k_terminal = fenics_convert(loadmat(f'{root_path}/../../data/terminal_state.mat', squeeze_me=True)["k_terminal"],"function", fun_space=K)

# Mass / stiffness matrices (block‑diagonal in time)
M_one = assemble(TrialFunction(K) * TestFunction(K) * dx).array()
M = block_diag(*[M_one] * num_steps)
K_mat_one = assemble(TrialFunction(K).dx(0) * TestFunction(K).dx(0) * dx).array()
K_mat = block_diag(*[K_mat_one] * num_steps)


# -------------------------------------------------------------------------
# Helper utilities
# -------------------------------------------------------------------------

def _manage_annotation(annotate: bool, verbose: bool) -> Tuple[int, Callable]:
    """
    Prepare the AD-tape and logging level.
    Returns the previous log level and a ``restore`` callable.
    """
    if annotate: get_working_tape().clear_tape()
    else: pause_annotation()

    prev_level = get_log_level()
    set_log_active(verbose)

    def restore():
        if not annotate: continue_annotation()
        set_log_level(prev_level)

    return prev_level, restore


def _time_step_loop(k0: Function,callback: Callable[[int, Function, Function], None],annotate: bool,verbose: bool) -> None:
    """
    Generic time-stepping loop.
    - ``k0`` : initial tracer state (Function in space K)
    - ``callback`` : called each step with (step_index, current_k, u_n)
    """
    _, restore = _manage_annotation(annotate, verbose)

    k = Function(K)          # current solution
    v_k = TestFunction(K)    # test function
    k_n = Function(K)
    k_n.assign(k0)

    t.assign(0.0)
    for n in range(num_steps):
        t.assign(float(t) + float(dt))

        # retrieve velocity at current time
        u_n = Function(U)
        u_timeseries.retrieve(u_n.vector(), float(t))

        # variational form (identical to original)
        F_k = (1 / dt * (k - k_n) * v_k + gamma * k.dx(0) * v_k.dx(0) + u_n * k.dx(0) * v_k + u_n.dx(0) * k * v_k + reac_fn(k) * v_k) * dx - k.dx(0) * v_k * ds
        solve(F_k == 0, k, J=derivative(F_k, k))
        k_n.assign(k)

        # user-defined processing for this step
        callback(n, k_n, u_n)

    restore()

# -------------------------------------------------------------------------
# Core model routines
# -------------------------------------------------------------------------

def state_solve(k0_input,return_type: Literal["vertex", "vector", "petsc", "function"],plot_k=False,annotate=True,verbose=False,return_all=False):
    """Forward solve for a single initial tracer field."""
    k0 = fenics_convert(k0_input, "function", K)

    # container for all intermediate states (if requested)
    states = []

    def step_callback(_, k_n, __):
        states.append(k_n.vector()[:])
        if plot_k: plot(k_n)

    _time_step_loop(k0, step_callback, annotate, verbose)

    if return_all: return np.array(states).flatten()
    return fenics_convert(states[-1], return_type, K)


def state_solve_all_obj(k0, kt_in, annotate=True, verbose=False):
    """Objective value for a given initial tracer and terminal data."""
    _, restore = _manage_annotation(annotate, verbose)

    # Prepare data
    k0_f = fenics_convert(k0, "function", K)
    k_n = Function(K)
    k_n.assign(k0_f)

    J_output = 0.0
    kt_in_mat = np.linalg.solve(M_one, fenics_convert(kt_in, "vector", K).reshape(num_steps, N+1).T)

    def step_callback(i, k_n, __):
        nonlocal J_output
        J_output += assemble(inner(k_n, fenics_convert(kt_in_mat[:, i], "function", K)) * dx)
        # inner product with the corresponding terminal test vector

    _time_step_loop(k0_f, step_callback, annotate, verbose)

    restore()
    return J_output


def J(k0, kt):
    """Regularised misfit functional."""
    k0_f = fenics_convert(k0, "function", K)
    kt_f = fenics_convert(kt, "function", K)
    misfit = 0.5 * inner(kt_f - k_terminal, kt_f - k_terminal) * dx
    reg = 0.5 * beta * inner(k0_f.dx(0), k0_f.dx(0)) * dx
    return 1e4 * assemble(misfit + reg)


def Jz(k0, kt):
    """Gradient of J with respect to z"""
    return 1e4 * float(beta) * K_mat_one @ fenics_convert(k0, "vector", K)


# -------------------------------------------------------------------------
# Reduced functional for optimization
# -------------------------------------------------------------------------

J_hat_np = None


def reduced_functional_J_hat(k0):
    """Create a NumPy-compatible reduced functional (cached)."""
    global J_hat_np
    k0_f = fenics_convert(k0, "function", K)
    kt = state_solve(k0_f, "function")
    J_val = J(k0_f, kt)
    control = Control(k0_f)
    J_hat = ReducedFunctional(J_val, control)
    J_hat_np = ReducedFunctionalNumPy(J_hat)
    return J_hat_np


def callback_call(J_hat, full=True):
    """Simple iteration callback for optimization algorithms."""
    c_test = Function(K)

    def callback_fn(x):
        if callback_fn.iteration == 0:
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


# -------------------------------------------------------------------------
# Interfaces used by the Matlab wrapper
# -------------------------------------------------------------------------

def apply_solution_operator_z_jacobian_transpose(kt_in, k0):
    """Adjoint of the forward operator (Jacobian^T) applied to a data vector."""
    k0_f = fenics_convert(k0, "function", K)
    grad = compute_gradient(state_solve_all_obj(k0_f, kt_in), Control(k0_f))
    return grad.vector()[:]


def misfit_gradient(kt, k0):
    """Gradient of the data-misfit term (compatible with padded vectors)."""
    kt_vec = fenics_convert(kt, "vector", K) 
    if kt_vec.size == N + 1: return 1e4 * M_one @ (kt_vec - fenics_convert(k_terminal, "vector", K))
    # padded case
    unpadded = 1e4 * M_one @ (kt_vec[-N - 1:] - fenics_convert(k_terminal, "vector", K))
    return np.pad(unpadded, (kt_vec.size - N - 1, 0), "constant")


def apply_misfit_hessian(kt_in, k0, kt):
    """Action of the misfit Hessian (constant, thus mass-scaled)."""
    kt_in_vec = fenics_convert(kt_in, "vector", K)
    if kt_in_vec.size == N + 1: raise Exception("Incorrect kt_in size!")
    # padded case
    unpadded = 1e4 * M_one @ fenics_convert(kt_in_vec[-N - 1:], "vector", K)
    return np.pad(unpadded, (kt_in_vec.size - N - 1, 0), "constant")

def apply_rs_hessian(k0_in, k0):
    """Reduced-space Hessian action for multiple right-hand sides."""
    k0_vec = fenics_convert(k0, "vector", K)
    k0_in_vec = fenics_convert(k0_in, "vector", K)
    if k0_in_vec.ndim == 1: k0_in_vec = k0_in_vec[:, np.newaxis]

    # Build the reduced functional lazily if not yet available
    if J_hat_np is None: reduced_functional_J_hat(k0_vec)
    return np.column_stack([J_hat_np.hessian(k0_vec, k0_in_vec[:, col])for col in range(k0_in_vec.shape[1])])

def apply_solution_operator_z_jacobian(k0_in, k0):
    """Jacobian of the forward operator applied to a perturbation."""
    return state_solve_all_jac(k0, k0_in)

def state_solve_all_jac(k0, k0_in, annotate=True, verbose=False):
    """Jacobian action for a single perturbation (used by OED)."""
    k0_f = fenics_convert(k0, "function", K)
    k0_in_f = fenics_convert(k0_in, "function", K)
    J_list = []

    def step_callback(_, k_n, __):
        J_list.append(compute_jacobian_action(k_n, Control(k0_f), k0_in_f).vector()[:])
    _time_step_loop(k0_f, step_callback, annotate, verbose)
    return np.array(J_list).flatten()


print(f"Succesfully imported {__name__}")