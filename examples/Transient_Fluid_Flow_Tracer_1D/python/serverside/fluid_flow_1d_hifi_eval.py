# IF PKG_CONFIG NOT FOUND
# try:
#     import matlab_env_fix
#     import os
#     os.environ = matlab_env_fix.data
#     import sys
#     sys.setdlopenflags(10)
# except ModuleNotFoundError:
#     pass
from fenics_helpers_noadj import *
from pathlib import Path
root_path = Path(__file__).parent
u_timeseries = TimeSeries(f"{root_path}/../../data/velocity_timeseries_midfi_1d")

N = 30;
unit_mesh = UnitIntervalMesh(N);
ds = generate_ds(unit_mesh)

# Define Function Spaces
P1 = FiniteElement('CG', unit_mesh.ufl_cell(), 1)
P2 = FiniteElement('CG', unit_mesh.ufl_cell(), 2)

C = FunctionSpace(unit_mesh, P1) # Tracer
U = FunctionSpace(unit_mesh, P2) # Velocity
P = FunctionSpace(unit_mesh, P1) # Pressure
R = FunctionSpace(unit_mesh, P2) # Density
E = FunctionSpace(unit_mesh, P2) # Internal Energy

# Mixed Function Space
CUPRE = FunctionSpace(unit_mesh, MixedElement([C.ufl_element(), U.ufl_element(), P.ufl_element(), R.ufl_element(), E.ufl_element()]))

# Obtain trial and test functions from mixed function space
cupre = Function(CUPRE)
c, u, p, r, e = split(cupre)
v_c, v_u, v_p, v_r, v_e = TestFunctions(CUPRE)
t = Constant(0);

# Set boundary conditions for u(x = 0) and u(x = 1)
# bc_c_left = DirichletBC(CUPRE.sub(0), Expression("exp(-(t+1))/pow(t+1, 2)", t=t, degree=2), "near(x[0], 0)")
# bc_c_right = DirichletBC(CUPRE.sub(0), Expression("2*exp(-(t+1))/pow(t+1, 2)", t=t, degree=2), "near(x[0], 1)")

# Set boundary conditions for u(x = 0) and u(x = 1)
bc_u_left = DirichletBC(CUPRE.sub(1), Expression("2/(t+1)", t=t, degree=2), "near(x[0], 0)")
bc_u_right = DirichletBC(CUPRE.sub(1), Expression("1/(t+1)", t=t, degree=2), "near(x[0], 1)")

# Set boundary conditions for p(x = 0) and p(x = 1)
# bc_p_left = DirichletBC(CUPRE.sub(2), Expression("1/pow(t+1, 2)", t=t, degree=2), "near(x[0], 0)")
# bc_p_right = DirichletBC(CUPRE.sub(2), Expression("2/pow(t+1, 2)", t=t, degree=2), "near(x[0], 1)")

# Set boundary conditions for r(x = 0) and r(x = 1)
bc_r_left = DirichletBC(CUPRE.sub(3), Constant(1), "near(x[0], 0)")
bc_r_right = DirichletBC(CUPRE.sub(3), Constant(1), "near(x[0], 1)")

# Set boundary conditions for e(x = 0) and e(x = 1)
# bc_e_left = DirichletBC(CUPRE.sub(4), Expression("1/pow(t+1, 2)", t=t, degree=2), "near(x[0], 0)")
# bc_e_right = DirichletBC(CUPRE.sub(4), Expression("2/pow(t+1, 2)", t=t, degree=2), "near(x[0], 1)")

# Combine Boundary Conditions
bcs = [bc_u_left, bc_u_right, bc_r_left, bc_r_right];

# Set Initial Condition at t = 0

rcv = Constant(1)
c0_exp = Expression("0", degree=2)
u0_exp = Expression("2-x[0]", degree=2)
r0_exp = Constant(1)
e0_exp = Expression("1*(1+x[0])", degree=2)
p0_exp = Expression("rcv*(1+x[0])", rcv=float(rcv), degree=2)

# u0_exp = Expression("x[0]+1", degree=2)
# p0_exp = Constant(1) # doesn't matter
# r0_exp = Expression("1/(1+x[0])", degree=2)
# e0_exp = Expression("1+9*x[0]", degree=2)

# Initialize u_n and p_n
c_n = interpolate(c0_exp, C)
u_n = interpolate(u0_exp, U)
p_n = interpolate(p0_exp, P) # just a newton guess
r_n = interpolate(r0_exp, R)
e_n = interpolate(e0_exp, E)

# Start from nonzero initial condition
ic = [c_n, u_n, p_n, r_n, e_n];
assign(cupre, ic)

# Set Grid and PDE Parameters
T = 0.1
num_steps = 25
dt = Constant(T/num_steps)
gamma = Constant(0.05)
# reac_fn = lambda c: Constant(2) * (c+Constant(1))**2
reac_fn = lambda c: Constant(1) * c
reac_enthalpy = Constant(100_000) # NOTE: MODIFIED!!!
kcv = Constant(30)

# Weak Form of PDE 
# F_3 = (1/dt*(e-e_n)*r*v_e + u*e.dx(0)*r*v_e + p*u.dx(0)*v_e + alpha*e.dx(0)*(r*v_e.dx(0) + r.dx(0)*v_e) ) * dx
# F_c = (1/dt * (c - c_n) * v_c + gamma*c.dx(0)*v_c.dx(0) + u*c.dx(0)*v_c + u.dx(0)*c*v_c + reac_fn(c)*v_c) * dx 
F_1 = (1/dt*(r-r_n)*v_r + u*r.dx(0)*v_r + r*u.dx(0)*v_r) * dx;
F_2 = (1/dt*(u-u_n)*r*v_u + u*u.dx(0)*r*v_u - p*v_u.dx(0)) * dx;
F_3 = (1/dt*(e-e_n)*r*v_e + u*e.dx(0)*r*v_e + p*u.dx(0)*v_e + kcv*e.dx(0)*v_e.dx(0) - reac_enthalpy*reac_fn(c)*v_e) * dx
F_4 = ((p - rcv*r*e)*v_p) * dx;
F_c = (1/dt * (c - c_n) * v_c + gamma*c.dx(0)*v_c.dx(0) + u*c.dx(0)*v_c + u.dx(0)*c*v_c + reac_fn(c)*v_c) * dx #- gamma*c.dx(0)*v_c*ds (imposing c.dx(0) = 0)
F = F_1 + F_2 + F_3 + F_4 + F_c

assemble(F); # This is just to get FFC JIT started well before PDE solve.

def state_solve(c0, return_type: Literal["vertex", "vector", "petsc", "function"], plot_c=False, return_all=False, store_midfi=False):
    global t

    # Reset initial conditions
    k_list = [];
    c_n.assign(fenics_convert(c0, "function", C))
    u_n = interpolate(u0_exp, U)
    p_n = interpolate(p0_exp, P) # just a newton guess
    r_n = interpolate(r0_exp, R)
    e_n = interpolate(e0_exp, E)
    ic = [c_n, u_n, p_n, r_n, e_n];
    assign(cupre, ic)
    if plot_c: plot(c_n)

    if store_midfi and len(u_timeseries.vector_times()) != 0: 
        store_midfi = False
        print("Cannot overwrite existing file.")

    u_temp = Function(U)

    # Solve the PDE with time-stepping
    t.assign(0)
    for n in range(num_steps):
        t.assign(float(t)+float(dt))
        
        # Solve Nonlinear PDE
        J = derivative(F, cupre)
        problem = NonlinearVariationalProblem(F, cupre, bcs, J)
        solver = NonlinearVariationalSolver(problem)
        solver.solve()

        # Assign cupre to concurrent iterates
        c_n.assign(cupre.sub(0, True))
        u_n.assign(cupre.sub(1, True))
        p_n.assign(cupre.sub(2, True))
        r_n.assign(cupre.sub(3, True))
        e_n.assign(cupre.sub(4, True))
        if store_midfi: u_timeseries.store(u_n.vector(), float(t))
        if plot_c: plot(c_n)
        if return_all: k_list.append(c_n.vector()[:])


        # u_timeseries.retrieve(u_temp.vector(), float(t))
        # print(np.linalg.norm(u_n.vector()[:]-u_temp.vector()[:])/np.linalg.norm(u_n.vector()[:]))
        # print(fenics_convert(u_n, 'vertex', U))
        # print(fenics_convert(u_temp, 'vertex', U))
        # print(np.linalg.norm(fenics_convert(u_n, 'vertex', U)))
        # print(np.linalg.norm(u_temp.vector()[:]))
    # print("Concentration: ", (fenics_convert(c_n, 'vertex', C)))
    # print("Velocity: ", (fenics_convert(u_n, 'vertex', U)))
    # print("Pressure: ", (fenics_convert(p_n, 'vertex', P)))
    # print("Density: ", (fenics_convert(r_n, 'vertex', R)))
    # print("Energy: ", (fenics_convert(e_n, 'vertex', E)))

    if return_all: return np.array(k_list).flatten()
    return fenics_convert(c_n, return_type);


# vec_grad = zeros(31, 1)
# J_z_bar = Jhat_hifi_fn(z_bar)
# for i=1:31
#     disp(i)
#     vec_grad(i) = 1e4*(Jhat_hifi_fn(z_bar+1e-4*I(:, i)) - J_z_bar)
# end