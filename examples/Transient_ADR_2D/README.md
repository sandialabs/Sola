# Two-dimensional Transient Advection-Diffusion

This example looks at a control problem for a set of two-species advection-diffusion-reaction dynamics in a two-dimensional urban canyon.
The physical situation is a contaminant ($u_1$) which will be contained by deploying a neutralizing retardant ($u_2$).

## Dynamics

Denoting the spatial domain as $\Omega$, the dynamics are the following:

$$
\tag{1}
\begin{aligned}
    \frac{\partial u_1}{\partial t}
    &= \kappa_1\Delta u_1 - \alpha_2 \mathbf{v}\cdot \nabla u_1 - \rho_1 u_1 u_2,
    \\
    \frac{\partial u_2}{\partial t}
    &= \kappa_2\Delta u_2 - \alpha_2 \mathbf{v}\cdot \nabla u_2 - \rho_2 u_1 u_2 + f,
    \\
    \vec{n}\cdot\nabla u_1 &= \vec{n}\cdot\nabla u_2 = 0
    ~~\textrm{ for }~~
    \mathbf{x}\in\partial\Omega,
    \\
    u_1(\mathbf{x}, 0) &= u_0(\mathbf{x}),
    \quad
    u_2(\mathbf{x}, 0) = 0,
\end{aligned}
$$

where $\mathbf{v}\in\mathbb{R}^{2}$ is a **constant** velocity field, $\kappa_1,\kappa_2 > 0$ are the diffusion coefficients, $\alpha_1,\alpha_2 > 0$ are the advection coefficients, $\rho_1,\rho_2 > 0$ are the reaction coefficients, and $f = f(\mathbf{x}, t)$ is a source term given by a sum of Gaussians:

$$
\begin{aligned}
    f(\mathbf{x},t)
    = \sum_{i=1}^{n_q}q_{i}(t)\phi_i(\mathbf{x}),
    \qquad\qquad
    \phi_i(\mathbf{x}) = 50 \exp\left(-1000\|\mathbf{x} - \bar{\mathbf{x}}_i\|^2\right),
\end{aligned}
$$

in which $\bar{\mathbf{x}}_i\in\mathbb{R}^2$ denotes the spatial center of the $i$-th Gaussian.

The (time-discrete) controls in this problem are the values of the source nodes $q_1(t),\ldots,q_m(t)$ at the $n_t - 1$ points of the temporal discretization, excluding the initial time point (because the eventual model is solved with the first-order implicit Euler time stepping scheme).

## Control Problem

For a mesh with $n_x\in\mathbb{N}$ nodes, let $\mathbf{y}_1(t),\mathbf{y}_2(t)\in\mathbb{R}^{n_x}$ denote the spatial discretizations of the PDE states $u_1$ and $u_2$, respectively, at time $t$.
The full discrete state is $\mathbf{y}(t) = [~\mathbf{y}_1(t)^{\mathsf{T}}~~\mathbf{y}_{2}(t)^{\mathsf{T}}~]^{\mathsf{T}}\in\mathbb{R}^{n_y}$ with $n_y = 2n_x$ degrees of freedom.
Let $\mathbf{q}(t)\in\mathbb{R}^{n_q}$ collect the source node coefficients at time $t$.
We consider the following minimization:

$$
\begin{aligned}
    &\min_{\mathbf{y}(t), \mathbf{q}(t)}
    \frac{1}{2}\int_{0}^{T}\|\mathbf{y}_1(t) \ast \mathbf{p}\|_{\mathbf{M}}^2\,dt
    + \frac{\gamma}{2}\int_{t_2}^{T}\|\mathbf{q}(t)\|_2^2\,dt
    \\
    &\textrm{subject to }~~(\mathbf{y}(t),\mathbf{q}(t))~~\textrm{ jointly solves the PDE (1)}
    \\
    &\textrm{and }~~q_i(t) \ge 0~~\textrm{ for all }~~i=1,\ldots,m~~\textrm{ and }~~t > 0.
\end{aligned}
$$

Here, $\mathbf{p}\in\mathbb{R}^{n_x}$ represents weights over the spatial domain indicating which areas should be prioritized to be protected from the contaminant, $\mathbf{M}\in\mathbb{R}^{n_x\times n_x}$ is a mass matrix, $\ast$ denotes the Hadamard (elementwise) product, and $\gamma > 0$ is a user-specified regularization constant.
The nonnegativity bound on the entries of $\mathbf{q}(t)$ ensures that $f$ is a source (adding retardant in), never a sink (taking retardant out).

The [`Transient_ADR_2D_Objective`](./Transient_ADR_2D_Objective.m) class implements this objective, discretizing the time integrals with the trapezoidal rule as usual and setting $\mathbf{p}$ as a Gaussian.
Note that the objective function only depends on $\mathbf{y}_1(t)$ (the contaminant) and $\mathbf{u}(t)$ (the input), not on the retardant.
We could penalize the amount of retardant as well (to model only having a limited amount of retardant), but a term like

$$
\begin{aligned}
    \int_{0}^{T}\|\mathbf{y}_2(t)\|_{\mathbf{M}}^2\,dt
\end{aligned}
$$

would penalize the total amount of retardant (which decreases as it hits the contaminant), not the amount pumped into the system via the input.

Alternatively, we could penalize the controller with the term

$$
\begin{aligned}
    \frac{\gamma}{2}\int_{t_2}^{T}\|\mathbf{Bq}(t)\|_\mathbf{M}^2\,dt
\end{aligned}
$$

where $\mathbf{B} = [~\phi_1(\mathbf{x})~~\cdots~~\phi_m(\mathbf{x})~] \in\mathbb{R}^{n_x \times n_q}$.
This would be perhaps a better representation of the actual source term, but it requires more intrusive information (the $\mathbf{B}$ matrix) and has the same overall effect of penalizing the entries of $\mathbf{q}(t)$.

Let $g(\mathbf{y}) = \|\mathbf{y}_1(t) \ast \mathbf{p}\|_{\mathbf{M}}^2 = (\mathbf{y}_1(t) \ast \mathbf{p})^{\mathsf{T}}\mathbf{M}(\mathbf{y}_1(t) \ast \mathbf{p}).$
Assuming $\mathbf{M}$ is symmetric, these are the derivatives of $g$:

$$
\begin{aligned}
    \nabla_{\!\mathbf{y}}g(\mathbf{y})
    = \left[\begin{array}{c}
        \nabla_{\!\mathbf{y}_1}g(\mathbf{y}) \\ \nabla_{\!\mathbf{y}_2}g(\mathbf{y})
    \end{array}\right]
    = \left[\begin{array}{c}
        \mathbf{p}\ast(\mathbf{M}(\mathbf{y}_1\ast\mathbf{p}))
        \\ \mathbf{0}
    \end{array}\right],
    \\
    \nabla_{\!\!\mathbf{y}}g(\mathbf{y})
    = \left[\begin{array}{cc}
        \nabla_{\!\!\mathbf{y}_1}g(\mathbf{y})
        & \mathbf{0} \\ \mathbf{0} & \mathbf{0}
    \end{array}\right]
    = \left[\begin{array}{cc}
        \mathbf{p}^{\mathsf{T}}\,\hat{\ast}\,\mathbf{M}\,\hat{\ast}\,\mathbf{p}
        & \mathbf{0} \\ \mathbf{0} & \mathbf{0}
    \end{array}\right],
\end{aligned}
$$

where $\hat{\ast}$ denotes broadcasted multiplication.

## Custom Mesh Generation

We've defined a rudimentary urban canyon mesh (two-dimensional with some holes) using MATLAB's PDE Modeler tool.
This section describes how to modify, export, and use the mesh.

- Within MATLAB, start the PDE Modeler with `> pdeModeler`
- Load the existing mesh: `File > Open...` and select the file (e.g., `canyon.m`).
- Make edits: add new boundaries, specify boundary conditions, set the mesh size, and so on.
- Export the boundary geometry and conditions: `Boundary > Export Decomposed Geometry, Boundary Cond's...`, then type `geometry bcs` in the dialogue box. This adds variables called `geometry` and `bcs` to the current MATLAB workspace.
- Export the mesh: `Mesh > Export Mesh...` and type `points edges triangles` in the dialogue box. This adds `points`, `edges`, and `triangles` to the current MATLAB workspace.
- Save the exported variables (below).

```matlab
> save('meshfile.mat', 'geometry', 'bcs', 'points', 'edges', 'triangles');
```

The static method `Transient_ADR_2D.model_from_file()` loads data from a `.mat` file and constructs a PDE model using the PDE Toolkit (which is different than the PDE Modeler in some ways).
Once the model is initialized, use `pdegplot(model,EdgeLabels="on")` to show the geometry and edge labels.

## Full-order Solver

The [`Transient_ADR_2D`](./Transient_ADR_2D.m) class uses MATLAB's [pde toolbox](https://www.mathworks.com/help/pde/ug/equations-you-can-solve.html) to solve the PDE (1) on the two-dimensional finite element mesh generated in the previous section.
See the script [Driver_Data_Generation.m](./Driver_Data_Generation.m) for examples of generating and visualizing full-order solutions.

We can extract a mass matrix from this solver, but not the differential operators due to the velocity, reaction, and source terms.

To animate the solution instead of redoing the computation:

```matlab
load('solver.mat', 'solver');
load('solution.mat', 'u');
solver.Animate_Solution(u.NodalSolution);
```

## Reduced-order Modeling

The script [Driver_Opt_OpInf_Multi.m](./Driver_Opt_OpInf_Multi.m) uses the `Transient_ADR_2D` solver to generate training data, learn a [POD basis](../../src/model_reduction/POD_Basis.m), and calibrate a reduced-order model through [Operator Inference](../../src/model_reduction/OpInf_ROM_Constraint_Multi.m):

$$
\begin{aligned}
    \frac{\textrm{d}}{\textrm{d}t}\hat{\mathbf{y}}_1(t)
    &= \hat{\mathbf{A}}_1\hat{\mathbf{y}}_1(t) + \hat{\mathbf{H}}_1[\hat{\mathbf{y}}_1\otimes\hat{\mathbf{y}}_2],
    \\
    \frac{\textrm{d}}{\textrm{d}t}\hat{\mathbf{y}}_2(t)
    &= \hat{\mathbf{A}}_2\hat{\mathbf{y}}_2(t) + \hat{\mathbf{H}}_2[\hat{\mathbf{y}}_1\otimes\hat{\mathbf{y}}_2] +
    \hat{\mathbf{B}}[\mathbf{q}(t) \ast \mathbf{q}(t)],
    \\
    \mathbf{y}_i(t)
    &\approx \mathbf{V}_{\!i}\hat{\mathbf{y}}_{i}(t), \quad i=1,2.
\end{aligned}
$$

Here, $\mathbf{V}_{i}\in\mathbb{R}^{n_x \times r_i}$ is a rank-$r_i$ POD basis learned from $\mathbf{y}_i(t)$ data.
Note that $r_1$ and $r_2$ do not have to be equal.
We then have reduced states $\mathbf{y}_{i}(t) \in \mathbb{R}^{r_i}$ and "operators" $\hat{\mathbf{A}}_i\in\mathbb{R}^{r_i\times r_i}$ and $\hat{\mathbf{H}}_i\in\mathbb{R}^{r_i \times r_1 r_2}$ for $i=1,2$, as well as $\hat{\mathbf{B}}\in\mathbb{R}^{r_2 \times n_q}$.
The Hadamard product $\ast$ is employed to ensure positivity in the control.

The regression is populated with sixth-order finite difference estimates of the time derivatives of the training states.

**Remark**. We could specify $\hat{\mathbf{B}}$ intrusively by constructing the control matrix $\mathbf{B}$ and setting $\hat{\mathbf{B}} = \mathbf{V}_{2}^\mathsf{T}\mathbf{B}$.
Likewise, since the nonlinearity $u_1 u_2$ is local, we can construct each $\hat{\mathbf{H}}_i$ intrusively if we like.
We want $\hat{\mathbf{H}}_i$ such that

$$
\begin{aligned}
    \hat{\mathbf{H}}_{i}[\hat{\mathbf{y}}_1\otimes\hat{\mathbf{y}}_2]
    &= \mathbf{V}_{i}^\mathsf{T}((\mathbf{V}_{\!1}\hat{\mathbf{y}}_1)\ast(\mathbf{V}_{\!2}\hat{\mathbf{y}}_2))
\end{aligned}
$$

It turns out that

$$
\begin{aligned}
    (\mathbf{V}_{\!1}\hat{\mathbf{y}}_1)\ast(\mathbf{V}_{\!2}\hat{\mathbf{y}}_2)
    = (\mathbf{V}_{\!1}^{\mathsf{T}} \odot \mathbf{V}_{\!2}^{\mathsf{T}})^\mathsf{T}[\hat{\mathbf{y}}_1\otimes\hat{\mathbf{y}}_2],
\end{aligned}
$$

where $\odot$ is the Khatri-Rao product.
Hence, $\hat{\mathbf{H}}_i = \mathbf{V}_i^\mathsf{T}(\mathbf{V}_{\!1}^{\mathsf{T}} \odot \mathbf{V}_{\!2}^{\mathsf{T}})^\mathsf{T}$.

(This is a mixed-product property, proven in Theorem 1 of [this paper](https://link.springer.com/article/10.1007/BF02733426)).

## Detailed Guide of `Driver_Opt_OpInf_Multi.m`

This section explains [Driver_Opt_OpInf_Multi.m](./Driver_Opt_OpInf_Multi.m), a script that

1. Generates full-order solutions to be used as training data,
2. Trains an OpInf ROM with the `Transient_ADR_2D_OpInf_Constraint` class,
3. Solves an optimization problem using the ROM as a surrogate, and
4. Visualizes results.

### Preliminaries

```matlab
clear;
close all;
clc;

%% Experiment parameters.

meshfile = 'urban_canyon.mat';
datafile = 'OpInf_Training_Data.mat';

regenerate_data = false;
plot_basis_functions = false;
plot_training_data = false;
plot_training_reconstruction = false;

residual_energies = [1e-3];
ABregularization_candidates = [0, 1, 10, 100];
Hregularization_candidates = logspace(2, 6, 21);
ddt_strategy = '6thOrder';
control_regularization = 5e-2;
```

Data variables:

- `meshfile`: file that contains the 2D geometry of the problem, [generated earlier](#custom-mesh-generation) with `pdeModeler`.
- `datafile`: file to save the high-fidelity training data (state snapshots and control profile) to.

Script directives:

- `regenerate_data`: if `true`, run the high-fidelity solver even if the `datafile` already exists and overwrite the `datafile` with the new training data.
- `plot_basis_functions`: if `true`, visualize the first $\min\{r_1,r_2\}$ POD basis vectors over the 2D spatial domain. Each pair of basis vectors is displayed in a separate figure, so this will make lots of figures.
- `plot_training_data`: if `true`, animate the high-fidelity solves when each is produced and, later, plot the compressed data in the coordinates of the POD bases.
- `plot_training_reconstruction`: if `true`, animate the training states over the 2D domain after projection to the space spanned by the POD basis functions.

OpInf parameters:

- `residual_energies`: the number of POD basis vectors is based on the residual energy level (a singular value criterion). For each entry of this vector, an OpInf ROM is learned (including regularization selection) and the training error is reported. This is more of a verification step than anything: only the final ROM is used for the optimization.
- `ABregularization_candidates`: potential scalars to regularize the linear terms in the model.
- `Hregularization_candidates`: potential scalars to regularize the quadratic terms in the model.
- `ddt_strategy`: how to estimate the time derivatives of the training states for the operator inference.

Optimization parameters:

- `control_regularization`: the scalar $\gamma$ used to penalize the controls in the objective function.

### Generate full-order solutions

```matlab
%% Generate training data if needed.

if ~exist(datafile, 'file') || regenerate_data
    disp('Generating training data');

    tic();
    % Initial condition parameters.
    init_center = [.05; .85];
    num_solves = 5;

    % Input function parameters.
    control_nodes = [0.1 0.5
                     0.1 0.9
                     0.1 1.1
                     0.3 0.7
                     0.3 0.9
                     0.3 1.1
                     0.5 0.3
                     0.5 0.5
                     0.5 0.7
                     0.7 0.7
                     0.9 0.3
                     0.9 1.1
                     1.1 0.7
                     1.1 0.9]';
    n_q = size(control_nodes, 2);
```

- The initial conditions are a Gaussian blob centered in space at `init_center`. Each run uses the same initial conditions, but different random controls. See the `Initial_Condition()` and `Initial_Contaminant()` methods of the `Transient_ADR_2D` class.
- `num_solves` is the number of high-fidelity solves to use for training data.
- The controls are Gaussian blob sources centered in space at the coordinates given in `control_nodes`. Physically, these nodes mark the locations where we can deploy species 2 (the decontaminant). See `Transient_ADR_2D.SourceTerm()`.
- `n_q` is the number of control nodes.

```matlab
    % Time domain.
    t = linspace(0, .4, 101);
    n_t = length(t);
    n_z = (n_t - 1) * n_q;

    % Load spatial geometry and mesh.
    model = Transient_ADR_2D.model_fromfile(meshfile);
    n_x = size(model.Mesh.Nodes, 2);
    n_y = 2 * n_x;
    n_u = n_y * n_t;

    % Model and input parameters.
    diffusion = [0.10, 0.10];
    advection = [4.00, 4.00];
    reaction = 2;
    num_randcontrol_nodes = 4;
    randcontrol_nodes = linspace(t(1), t(end), num_randcontrol_nodes);
```

- The controls are measured at each time step except the initial time, so the total control dimension is $n_z = n_q (n_t - 1)$.
- There are two species, so the total state dimension is $n_y = 2n_x$ where $n_x$ is the number of spatial nodes.
- `diffusion`, `advection`, and `reaction` are the (fixed) scalar parameters in the governing equations. In the language of this document, `diffusion` is $(\kappa_1, \kappa_2)$, `advection` is $(\alpha_1,\alpha_2)$, and `reaction` is $\rho_1 = \rho_2$.
- The training trajectories have random control profiles constructed as random splines with `num_randcontrol_nodes` nodes. For each source term, each of the nodes will be assigned a random positive value.

```matlab
    Z_train = zeros(n_z, num_solves);
    U_train = zeros(n_u, num_solves);

    for k = 1:num_solves
        disp(['High-fidelity solve ', num2str(k)]);

        % Initialize the solver.
        solver = Transient_ADR_2D(model, init_center, ...
                                  diffusion, advection, reaction, control_nodes);

        % Set up a random control profile.
        vals = [zeros(n_q, 1), 50 * rand(n_q, num_randcontrol_nodes - 1)];
        pp = spline(randcontrol_nodes, vals);
        controller = @(tt) ppval(pp, tt);

        % Solve the system.
        Yk = solver.State_Solve(controller, t).NodalSolution;

        if plot_training_data
            solver.Animate_Solution(Yk);
        end

        Qk = controller(t(2:end));

        % Record results.
        U_train(:, k) = reshape(Yk, [], 1);
        Z_train(:, k) = reshape(Qk, [], 1);
    end
    time_trainingdata = toc();

    save(datafile, "t", "solver", "U_train", "Z_train", "time_trainingdata");
end
```

- `Yk` is the high-fidelity state solution and `Qk` is the corresponding control. Each are flattened and stored as columns of `U_train` and `Z_train`, respectively.

```matlab
%% Load training data.

load(datafile);
n_t = length(t);
T = t(end);
n_u = size(U_train, 1);
num_solves = size(U_train, 2);
n_y = n_u / n_t;
n_x = n_y / 2;  % = size(solver.model.Mesh.Nodes, 2);
mass_matrix = assembleFEMatrices(solver.model, 'M').M;
mass_matrix = mass_matrix(1:n_x, 1:n_x);
stiffness_matrix = assembleFEMatrices(solver.model, 'K').K;
stiffness_matrix = stiffness_matrix(1:n_x, 1:n_x);

n_z = size(Z_train, 1);
n_q = n_z / (n_t - 1);  % = solver.n_q;
fprintf('Using %d training trajectories\n', num_solves);
```

- This section loads all variables from the `datafile`, even if the file was just generated.
- The mass matrix resulting from `assembleFEMatrices()` is $n_y \times n_y$ (or $2n_x\times 2n_x$), but this is just one $n_x \times n_x$ matrix repeated twice in block diagonal form, so we pull out the top left block. Same for the stiffness matrix, which is not used in this script but which is needed for the model discrepancy.

```matlab
%% Learn a POD basis for each variable.

% Unpack the states and controls by training trajectory.
states = cell(num_solves);
controls = cell(num_solves);
controls_romtraining = cell(num_solves);
for k = 1:num_solves
    states{k} = reshape(U_train(:, k), n_y, n_t);
    controls{k} = reshape(Z_train(:, k), n_q, n_t - 1);
    controls_romtraining{k} = sqrt(abs(controls{k}));
end

% Learn POD bases from the collection of all state snapshots.
states_all = horzcat(states{:});
basis1 = POD_Basis(states_all(1:n_x, :), false);  % , full(mass_matrix));
basis1.Set_Reduced_Dimension_From_Residual_Energy(residual_energies(1));
basis2 = POD_Basis(states_all(n_x + 1:end, :), false);  % , full(mass_matrix));
basis2.Set_Reduced_Dimension_From_Residual_Energy(residual_energies(1));
```

- First we reshape the training trajectories from $n_u \times 1$ to $n_y \times n_t$ (call these $\mathbf{Y}_k$) and the controls from $n_z \times 1$ to $n_q \times (n_t - 1)$.
- `controls_romtraining` is the square root of the controls so that, if $\mathbf{q}(t)$ is `controls_romtraining`, $\mathbf{q}(t)\ast\mathbf{q}(t)$ recovers the original `controls`.
- Next we concatenate the training trajectories to $\mathbf{Y} = [~~\mathbf{Y}_1~\cdots~\mathbf{Y}_{n_\text{trajectories}}~~]$. We take SVD of the first $n_x$ rows of $\mathbf{Y}$ to form the POD basis for the first species, and the SVD of the last $n_x$ rows of $\mathbf{Y}$ for the POD basis for the second species.
- The mass matrix is neglected by default because inverting the mass matrix takes a long time, but it can be included by using `basis1 = POD_Basis(states_all(1:n_x, :), false, full(mass_matrix));` and similar for `basis2`.

```matlab
if plot_basis_functions
    for i = 1:min(basis1.r, basis2.r)
        solver.Plot_Field([basis1.V(:, i), basis2.V(:, i)]);
        title(['POD basis function ', num2str(i)]);
    end
end

if plot_training_data
    for k = 1:num_solves
        Yhatk_1 = basis1.Compress(states{k}(1:n_x, :));
        Yhatk_2 = basis2.Compress(states{k}(n_x + 1:end, :));
        Yhatk = [Yhatk_1; Yhatk_2];
        figure;
        plot(t, Yhatk);
        title(['compressed state training data, trajectory', num2str(k)]);
    end
end
```

- Basis functions are visualized over the 2D spatial domain.
- The training data are compressed in the POD space, then we plot the POD coefficients as a function of time.

### Train an OpInf ROM

The next block loops through the `residual_energies` and does OpInf for each choice of residual energy, which dictates the number of POD basis vectors for each state.

```matlab
%% Learn a ROM, varying the reduced state dimension.

errors = zeros(length(residual_energies), 1);
for i = 1:length(residual_energies)
    res_energy = residual_energies(i);
    fprintf('\nUsing %.2e residual energy\n', res_energy);

    basis1.Set_Reduced_Dimension_From_Residual_Energy(res_energy);
    basis2.Set_Reduced_Dimension_From_Residual_Energy(res_energy);
    r_1 = basis1.r;
    r_2 = basis2.r;
    n_yr = r_1 + r_2;
    fprintf('POD with r_1 = %d and r_2 = %d basis vectors\n', r_1, r_2);

    % Compress states and check projection error.
    states_lofi = cell(num_solves);
    for k = 1:num_solves
        Yhat_1 = basis1.Compress(states{k}(1:n_x, :));
        Yhat_2 = basis2.Compress(states{k}(n_x + 1:end, :));
        states_lofi{k} = [Yhat_1; Yhat_2];
        Yproj_1 = basis1.Decompress(Yhat_1);
        Yproj_2 = basis2.Decompress(Yhat_2);
        Yproj = [Yproj_1; Yproj_2];
        proj_err = norm(Yproj - states{k}) / norm(states{k});
        fprintf("Projection error of trajectory %d: %.4f%%\n", k, 100 * proj_err);
    end
```

- The states are compressed to the POD basis as `states_lofi`.
- As a check, the compressed states are decompressed and compared to the original states.

```matlab
    %% Learn an OpInf ROM from the data.

    rom = Transient_ADR_2D_OpInf_Constraint(r_1, r_2, n_q, T, n_t, zeros(n_yr, 1));
    tic();
    rom.Select_Regularization(states_lofi, controls_romtraining, ...
                              ABregularization_candidates, ...
                              Hregularization_candidates, ...
                              ddt_strategy);
    time_opinfcalibration = toc();
```

The `Select_Regularization` method estimates the time derivatives of the training states (`ddt_strategy`), uses each combination of entries from `ABregularization_candidates` and `Hregularization_candidates` to define Tikhonov regularizers, and solves the ROM for each set of training controls. The regularizer that produces the least training error is deemed the winner.

```matlab
    % Solve the ROM for each of the training controls.
    total_error = 0;
    for k = 1:num_solves
        Yk_data = states{k};
        rom.y0 = states_lofi{k}(:, 1);
        Yk_rom_compressed = rom.State_Solve2(controls_romtraining{k});
        Yk_rom_1 = basis1.Decompress(Yk_rom_compressed(1:r_1, :));
        Yk_rom_2 = basis2.Decompress(Yk_rom_compressed(r_1 + 1:end, :));
        Yk_rom = [Yk_rom_1; Yk_rom_2];
        state_error = norm(Yk_data - Yk_rom) / norm(Yk_data);
        fprintf('ROM reconstruction error for training set %d: %.2f%%\n', ...
                k, 100 * state_error);
        total_error = total_error + state_error;
        if plot_training_reconstruction
            solver.Animate_Solution(Yk_rom);
        end
    end
    errors(i) = total_error / num_solves;
end

if length(residual_energies) > 1
    figure;
    semilogx(residual_energies, errors);
    title('Residual energy versus average ROM training error');
end
```

- The best-regularization ROM is solved for each set of training controls and the results are compared to the true states. This is a check that the ROM is reasonable: if these numbers are bad, there is little hope for the optimization.
- At the end of the loop, we plot the average reconstruction error for each choice of residual energy (each basis size).
- Only the final ROM (the one with residual energy `residual_energies(end)`) is used in the next step for optimization.

### Optimization with the OpInf ROM

```matlab
%% Set up the optimization objective.

% Make sure the initial conditions are right.
solver.init_center = [.05; .85];
rom.y0 = states_lofi{1}(:, 1);

obj_hifi = solver.Make_Objective([.6; .6], t(end), length(t), control_regularization);
Vfull = blkdiag(basis1.V, basis2.V);
obj_lofi = Reduced_Dynamic_Objective(obj_hifi, Vfull);
solver.Plot_Field(obj_hifi.target_weight, 'Protection zone');
```

- The `Reduced_Dynamic_Objective` class (see docs) modifies the high-fidelity objective for POD-style reduced states.
- We plot the section of the 2D domain that we want to protect from the first species (the contaminant).
- In `solver.Make_Objective`, `control_regularization` sets the strength of the regularizer on the control, $\gamma$. Increasing this number should increasingly penalize the total amount of decontaminant used.

```matlab
%% Set up and solve the optimization problem (using last trained ROM).

opt = Reduced_Space_Optimization(obj_lofi, rom);
% opt.z_lb = zeros(n_z, 1);                   % Lower bounds for control.
% opt.z_ub = 25 * ones(n_z, 1);               % Upper bounds for control.
opt.max_cg_iter = 200;

tic();
[u_lofi, z_lofi] = opt.Optimize(rand(n_z, 1));
time_lofioptimization = toc();
```

- `opt.z_lb` and `opt.z_ub` set lower and upper bounds, respectively, on the controls. The controls should be positive in this setup.
- This step will take a while, even with a ROM surrogate, but it shouldn't take more than an hour or so.

### Visualize Results

```matlab
%% Visualize optimization results.

% Inspect the state solution.
u_lofi_reshape = reshape(u_lofi, n_yr, n_t);
Y_rom_1 = basis1.Decompress(u_lofi_reshape(1:r_1, :));
Y_rom_2 = basis2.Decompress(u_lofi_reshape(r_1 + 1:end, :));
Y_rom = [Y_rom_1; Y_rom_2];
solver.Animate_Solution(Y_rom);             % ROM state with ROM controller

% Inspect the control solution.
Q_rom = reshape(z_lofi, n_q, n_t - 1).^2;
figure;
plot(t(2:end), Q_rom);
title('Optimal controls (optimized with an OpInf ROM surrogate)');
```

- `Y_rom` is the state solution to the surrogate-driven optimization problem, mapped back to the original state space.
- `Q_rom` is the control solution. Note that we square the ROM controls to get FOM controls. We plot the controls pointwise, which tells us how much each source is being turned on as a function of time (how much contaminant is being deployed at a particular location).

```matlab
% Solve the high-fidelity model with the inferred controls.
disp('Final high-fidelity solve');
pp = spline(t(2:end), Q_rom);
controller = @(tt) ppval(pp, tt);
Y_hifi = solver.State_Solve(controller, t).NodalSolution;
solver.Animate_Solution(Y_hifi);            % FOM state with ROM controller

save('OptimizationSolution.mat', "solver", "Y_hifi", "Y_rom", "t", "Q_rom", "n_q");

%% Load and visualize results later.
% load('OptimizationSolution.mat', "solver", "Y_hifi", "Y_rom", "t", "Q_rom", "n_q");
% figure;
% plot(t(2:end), Q_rom);
% title('Optimal controls (optimized with an OpInf ROM surrogate)');
% solver.Animate_Solution(Y_rom);   % ROM state with ROM controller
% solver.Animate_Solution(Y_hifi);  % FOM state with ROM controller
```

This final block solves the high-fidelity system with the control profile produced by the surrogate-driven optimization.
The results are saved for later, and you can load the results using the last bit of commented code.
