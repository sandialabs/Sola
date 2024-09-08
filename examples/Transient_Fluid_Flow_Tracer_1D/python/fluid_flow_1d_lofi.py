from fenics_helpers import * 
from pathlib import Path
from pyadjoint.reduced_functional_numpy import ReducedFunctionalNumPy
from scipy.io import loadmat

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
dt = Constant(T/num_steps)
t = Constant(0);

# # Retreive velocity from Timeseries
u_timeseries = TimeSeries(f"{root_path}/../data/velocity_timeseries_lofi_1d")

# Weak form of PDE
gamma = Constant(0.025)
reac_fn = lambda c: Constant(10) * c

# Store Mass Matrix for Future
M = assemble(TrialFunction(K) * TestFunction(K) * dx).array()
K_mat = assemble(TrialFunction(K).dx(0) * TestFunction(K).dx(0) * dx).array()


# Initial setup for inverse problem
k_terminal = fenics_convert(loadmat(f'{root_path}/../data/terminal_state.mat', squeeze_me=True)["k_terminal"], "function", fun_space=K)
beta = Constant(1e-5)

def state_solve(k0_input, return_type: Literal["vertex", "vector", "petsc", "function"], plot_k=False, annotate=True, verbose=False):
    # Handle annotation and verbosity
    if annotate: get_working_tape().clear_tape()
    else: pause_annotation()
    log_level = get_log_level()
    set_log_active(verbose)

    # Set Trial and Test Functions
    k = Function(K)
    v_k = TestFunction(K)

    # Set initial conditions
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
        F_k = (1/dt * (k - k_n) * v_k + gamma*k.dx(0)*v_k.dx(0) + u_n*k.dx(0)*v_k + u_n.dx(0)*k*v_k + reac_fn(k)*v_k) * dx #- gamma*k.dx(0)*v_k*ds
        J = derivative(F_k, k)
        solve(F_k == 0, k, J=J)
        k_n.assign(k)
        if plot_k: plot(k_n)

    # Revert annotation and verbosity
    if not annotate: continue_annotation()
    set_log_level(log_level)

    # Return final state
    return fenics_convert(k_n, return_type, K)

def J(k0, kt):
    # Convert inputs to functions
    k0 = fenics_convert(k0, "function", fun_space=K)
    kt = fenics_convert(kt, "function", fun_space=K)
    val = assemble(0.5*inner(kt - k_terminal, kt - k_terminal)*dx + 0.5 * beta * inner(k0.dx(0), k0.dx(0)) * dx)
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
def apply_solution_operator_z_jacobian_transpose(kt_in, k0):
    # Multiply by inv(M) to compensate for assembly in L2
    kt_in = np.linalg.solve(M, fenics_convert(kt_in, "vector", fun_space=K))
    kt_in = fenics_convert(kt_in, "function", fun_space=K)
    k0 = fenics_convert(k0, "function", fun_space=K)
    kt = state_solve(k0, return_type="function")
    return compute_gradient(assemble(inner(kt, kt_in) * dx), Control(k0)).vector()[:]

def misfit_gradient(kt, k0):
    return M @ (fenics_convert(kt, "vector", K)  - fenics_convert(k_terminal, "vector"));

def apply_misfit_hessian(kt_in, k0, kt):
    return M @ fenics_convert(kt_in, "vector", K)

def apply_rs_hessian(k0_in, k0):
    k0 = fenics_convert(k0, "vector")
    k0_in = fenics_convert(k0_in, "vector")
    if J_hat_np is None: reduced_functional_J_hat(k0); 
    if k0_in.ndim == 1: k0_in = k0_in[:, np.newaxis]
    return np.column_stack([J_hat_np.hessian(k0, k0_in[:, col]) for col in range(k0_in.shape[1])])

# Extra (just-in-case)
def apply_solution_operator_z_jacobian(k0_in, k0):
    k0_in = fenics_convert(k0_in, "function", fun_space=K)
    k0 = fenics_convert(k0, "function", fun_space=K)
    kt = state_solve(k0, return_type="function")
    return compute_jacobian_action(kt, Control(k0), k0_in).vector()[:]
        

print(f"Succesfully imported {__name__}")

