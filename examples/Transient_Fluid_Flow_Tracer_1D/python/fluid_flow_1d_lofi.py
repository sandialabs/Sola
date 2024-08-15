from fenics_helpers import * 
from scipy.io import loadmat
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
# Define Temporal Mesh
T = 0.1
num_steps = 25
dt = Constant(T/num_steps)
t = Constant(0);

# # Retreive velocity from Timeseries
u_timeseries = TimeSeries("data/velocity_timeseries_midfi_1d")

# Weak form of PDE
gamma = Constant(0.025)
reac_fn = lambda c: Constant(1) * c

# Store Mass Matrix for Future
M = assemble(TrialFunction(K) * TestFunction(K) * dx).array()
K_mat = assemble(TrialFunction(K).dx(0) * TestFunction(K).dx(0) * dx).array()
jacobian = None; jac_k0 = None; jac_kt = None;

# Initial setup for inverse problem
k_terminal = fenics_convert(loadmat('data/terminal_state.mat', squeeze_me=True)["k_terminal"], "function", fun_space=K)
beta = Constant(1e-5)

def state_solve(k0_input, return_type: Literal["vertex", "vector", "petsc", "function"] = "vertex", plot_k=False, annotate=False, verbose=True):
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

def J(k0_guess, k_n, return_gradients=True):
    # Convert inputs to functions
    k0_guess = fenics_convert(k0_guess, "function", fun_space=K)
    k_n = fenics_convert(k_n, "function", fun_space=K)

    val = assemble(0.5*inner(k_n - k_terminal, k_n - k_terminal)*dx + 0.5 * beta * inner(k0_guess.dx(0), k0_guess.dx(0)) * dx)
    if not return_gradients: return val
    grad_u = M @ (k_n.vector()[:] - k_terminal.vector()[:]);
    grad_z = float(beta) * K_mat @ k0_guess.vector()[:];
    return [val, grad_z, grad_u]

def J_uu_apply(k0_guess, k_n, v):
    return M @ fenics_convert(v, "vector", K)
    
def J_zz_apply(k0_guess, k_n, v):
    return float(beta) * K_mat @ fenics_convert(v, "vector", K)
    
def setup_inverse_problem(k0_guess, k_n, k_terminal_input=None, beta_input=None):
    # Set global beta & k_terminal
    global beta, k_terminal, J_hat
    if beta_input is not None: beta.assign(beta_input)
    if beta_input is not None: print("Deprecation of beta-input via setup_inverse_problem.")
    k0_guess = fenics_convert(k0_guess, "function", fun_space=K)
    k_n = fenics_convert(k_n, "function", fun_space=K)

    J_inv = J(k0_guess, k_n, return_gradients=False)
    control = Control(k0_guess)
    J_hat = ReducedFunctional(J_inv, control)
    return J_hat

def J_hat_hessian(k0, k0_in):
    k0 = fenics_convert(k0, "function", fun_space=K)
    k0_in = fenics_convert(k0_in, "function", fun_space=K)
    J_hat(k0);
    return fenics_convert(J_hat.hessian(k0_in), "vector")


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

def eval_c(k0, kt, return_type="function", fun_space=K):
    k0 = fenics_convert(k0, "function", fun_space=K)
    kt = fenics_convert(kt, "function", fun_space=K)
    c = Function(K)
    c.assign(state_solve(k0, return_type="function", verbose=False, annotate=True) - kt)
    return fenics_convert(c, return_type)

def compute_jacobian(c, control):
    global jacobian
    print("Building Jacobian [this may take some time]...")
    # Initialize an empty list to store the Jacobian rows
    jacobian_rows = []
    
    # Loop over each degree of freedom
    for i in range(K.dim()):
        # Create a unit vector in the direction of the i-th degree of freedom
        unit_vector = Function(K)
        unit_vector.vector()[i] = 1.0
        
        # Compute the derivative of the i-th component of c with respect to the control
        dci_dk = compute_gradient(assemble(inner(c, unit_vector) * dx), control)
        jacobian_rows.append(dci_dk.vector()[:])
    
    jacobian = np.linalg.solve(M, np.array(jacobian_rows)).T;

def c_z_transpose_apply_built(kt_in, k0, kt):
    global jac_k0, jac_kt # cached jacobian inputs
    k0 = fenics_convert(k0, "function", fun_space=K)
    kt = fenics_convert(kt, "function", fun_space=K)
    if jacobian is None or not np.allclose(jac_k0, k0.vector()[:]) or not np.allclose(jac_kt, kt.vector()[:]):
        c = eval_c(k0, kt)
        compute_jacobian(c, Control(k0));
        jac_k0 = k0.vector()[:]; jac_kt = kt.vector()[:]
    return jacobian.T @ fenics_convert(kt_in, "vector")

def c_z_transpose_apply(kt_in, k0, kt):
    kt_in = fenics_convert(kt_in, "function", fun_space=K)
    k0 = fenics_convert(k0, "function", fun_space=K)
    kt = fenics_convert(kt, "function", fun_space=K)
    c = eval_c(k0, kt)
    return compute_gradient(assemble(inner(c, kt_in) * dx), Control(k0)).vector()[:]

def c_z_apply(k0_in, k0, kt):
    global jac_k0, jac_kt # cached jacobian inputs
    k0 = fenics_convert(k0, "function", fun_space=K)
    kt = fenics_convert(kt, "function", fun_space=K)
    if jacobian is None or not np.allclose(jac_k0, k0.vector()[:]) or not np.allclose(jac_kt, kt.vector()[:]):
        c = eval_c(k0, kt)
        compute_jacobian(c, Control(k0))
        jac_k0 = k0.vector()[:]; jac_kt = kt.vector()[:]
    return jacobian @ fenics_convert(k0_in, "vector")

def c_u_inv_apply(kt_in, k0, kt):
    kt_in = fenics_convert(kt_in, "vector", fun_space=K)
    return - kt_in

def c_u_inv_transpose_apply(kt_in, k0, kt):
    kt_in = fenics_convert(kt_in, "vector", fun_space=K)
    return - kt_in
    
print(f"Succesfully imported {__name__}")
