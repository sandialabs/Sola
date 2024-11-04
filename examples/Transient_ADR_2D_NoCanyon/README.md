# Two-dimensional Transient Advection-Diffusion without Obstacles

This example looks at a control problem for a set of two-species advection-diffusion-reaction dynamics in a convex two-dimensional domain.
The physical situation is a contaminant ($u_1$) which will be contained by deploying a neutralizing retardant ($u_2$).

## Dynamics

Denoting the spatial domain as $\Omega$, the dynamics are the following:

$$
\tag{1}
\begin{aligned}
    \frac{\partial u_1}{\partial t}
    &= \kappa_1\Delta u_1 - \alpha \mathbf{v}\cdot \nabla u_1 - \rho u_1 u_2,
    \\
    \frac{\partial u_2}{\partial t}
    &= \kappa_2\Delta u_2 - \alpha \mathbf{v}\cdot \nabla u_2 - \rho u_1 u_2 + f,
    \\
    \vec{n}\cdot\nabla u_1 &= \vec{n}\cdot\nabla u_2 = 0
    ~~\textrm{ for }~~
    \mathbf{x}\in\partial\Omega,
    \\
    u_1(\mathbf{x}, 0) &= u^{(0)}_1(\mathbf{x}),
    \quad
    u_2(\mathbf{x}, 0) = 0,
\end{aligned}
$$

where $\mathbf{v} = [0, 1]^\mathsf{T}\in\mathbb{R}^{2}$ is the constant velocity field, $\kappa_1,\kappa_2 > 0$ are diffusion coefficients, $\alpha > 0$ is an advection coefficient, $\rho > 0$ is a reaction coefficient, and $f = f(\mathbf{x}, t)$ is a source term given by a sum of Gaussians:

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
    &\textrm{where }~~(\mathbf{y}(t),\mathbf{q}(t))~~\textrm{ jointly solves the PDE (1)}
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

We've defined a two-dimensional mesh using MATLAB's PDE Modeler tool.
This section describes how to modify, export, and use the mesh.

- Within MATLAB, start the PDE Modeler with `> pdeModeler`
- Load the existing mesh: `File > Open...` and select the file (e.g., `mesh.m`).
- Make edits: add new boundaries, specify boundary conditions, set the mesh size, and so on.
- Export the boundary geometry and conditions: `Boundary > Export Decomposed Geometry, Boundary Cond's...`, then type `geometry bcs` in the dialogue box. This adds variables called `geometry` and `bcs` to the current MATLAB workspace.
- Export the mesh: `Mesh > Export Mesh...` and type `points edges triangles` in the dialogue box. This adds `points`, `edges`, and `triangles` to the current MATLAB workspace.
- Save the exported variables (below).

```matlab
> save('newmesh.mat', 'geometry', 'bcs', 'points', 'edges', 'triangles');
```

The static method `Transient_ADR_2D.model_from_file()` loads data from a `.mat` file and constructs a PDE model using the PDE Toolkit (which is different than the PDE Modeler in some ways).
Once the model is initialized, use `pdegplot(model, EdgeLabels='on')` to show the geometry and edge labels.

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

The script [Driver_OptOpInf.m](./Driver_OptOpInf.m) uses the `Transient_ADR_2D` solver to generate training data, learn a [POD basis](../../src/model_reduction/POD_Basis.m), and calibrate a reduced-order model through [Operator Inference](../../src/model_reduction/OpInf_ROM_Constraint_Multi.m):

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
