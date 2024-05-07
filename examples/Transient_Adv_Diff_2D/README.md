# Two-dimensional Transient Advection-Diffusion

This example looks at a control problem for advection-diffusion dynamics in a two-dimensional urban canyon.

## Dynamics

Denoting the spatial domain as $\Omega$, the dynamics are the following:

$$
\tag{1}
\begin{aligned}
    \frac{\partial u}{\partial t}
    = \kappa\frac{\partial^{2} u}{\partial^2 x} - \alpha \mathbf{v}\cdot \nabla u + f,
    \\
    \vec{n}\cdot\nabla u = 0 ~~\textrm{ for }~~ \mathbf{x}\in\partial\Omega,
    \\
    u(\mathbf{x}, 0) = u_0(\mathbf{x}),
\end{aligned}
$$

where $\mathbf{v}\in\mathbb{R}^{2}$ is a **constant** velocity field, $\kappa > 0$ is the diffusion coefficient, $\alpha > 0$ is the advection coefficient, and $f = f(\mathbf{x}, t)$ is a source term given by a sum of Gaussians:

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

Let $\mathbf{y}(t)\in\mathbb{R}^{n_y}$ be a spatial discretization of the PDE state $u$ at time $t$, and let $\mathbf{q}(t)\in\mathbb{R}^{n_q}$ collect the source node coefficients at time $t$.
We consider the following minimization:

$$
\begin{aligned}
    &\min_{\mathbf{y}(t), \mathbf{q}(t)}
    \frac{1}{2}\int_{0}^{T}\|\mathbf{y}(t) \odot \mathbf{p}\|_{\mathbf{M}}^2\,dt
    + \frac{\gamma}{2}\int_{t_2}^{T}\|\mathbf{q}(t)\|_2^2\,dt
    \\
    &\textrm{subject to }~~(\mathbf{y}(t),\mathbf{q}(t))~~\textrm{ jointly solves the PDE (1)}
    \\
    &\textrm{and }~~q_i(t) \le 0~~\textrm{ for all }~~i=1,\ldots,m~~\textrm{ and }~~t > 0.
\end{aligned}
$$

Here, $\mathbf{p}\in\mathbb{R}^{n_y}$ represents weights over the spatial domain indicating which areas should be prioritized to be protected from the contaminant, $\mathbf{M}\in\mathbb{R}^{n_y\times n_y}$ is a mass matrix, $\odot$ denotes the Hadamard (elementwise) product, and $\gamma > 0$ is a user-specified regularization constant.
The nonpositivity bound on the entries of $\mathbf{q}(t)$ ensures that $f$ is a sink (taking contaminant out), never a source (adding contaminant in).

The [`Transient_Adv_Diff_2D_Objective`](./Transient_Adv_Diff_2D_Objective.m) class implements this objective, discretizing the time integrals with the trapezoidal rule as usual and setting $\mathbf{p}$ as a Gaussian.

Note that we could penalize the controller with the term

$$
\begin{aligned}
    \frac{\gamma}{2}\int_{t_2}^{T}\|\mathbf{Bq}(t)\|_\mathbf{M}^2\,dt
\end{aligned}
$$

where $\mathbf{B} = [~\phi_1(\mathbf{x})~~\cdots~~\phi_m(\mathbf{x})~] \in\mathbb{R}^{n_y \times n_q}$.
This would be perhaps a better representation of the actual source term, but it requires more intrusive information (the $\mathbf{B}$ matrix) and has the same effect of penalizing the entries of $\mathbf{q}(t)$.

## Custom Mesh Generation

We've defined a rudimentary urban canyon mesh (two-dimensional with some holes) using MATLAB's PDE Modeler tool.
This section describes how to modify, export, and use the mesh.

- Within MATLAB, start the PDE Modeler with `> pdeModeler`
- Load the existing mesh: `File > Open...` and select the file (e.g., `canyon.m`).
- Make edits: add new boundaries, specify boundary conditions, set the mesh size, and so on.
- Export the boundary geometry and conditions: `Boundary > Export Decomposed Geometry, Boundary Cond's...`, then type `decgeometry bcs` in the dialogue box. This adds variables called `geometry` and `bcs` to the current MATLAB workspace.
- Export the mesh: `Mesh > Export Mesh...` and type `points edges triangles` in the dialogue box. This adds `points`, `edges`, and `triangles` to the current MATLAB workspace.
- Save the exported variables (below).

```matlab
> save('meshfile.mat', 'geometry', 'bcs', 'points', 'edges', 'triangles');
```

The static method `Transient_Adv_Diff_2D.model_from_file()` loads data from a `.mat` file and constructs a PDE model using the PDE Toolkit (which is different than the PDE Modeler in some ways).
Once the model is initialized, use `pdegplot(model,EdgeLabels="on")` to show the geometry and edge labels.

## Full-order Solver

The [`Transient_Adv_Diff_2D`](./Transient_Adv_Diff_2D.m) class uses MATLAB's [pde toolbox](https://www.mathworks.com/help/pde/ug/equations-you-can-solve.html) to solve the PDE (1) on the two-dimensional finite element mesh generated in the previous section.
See the script [Driver_Data_Generation.m](./Driver_Data_Generation.m) for examples of generating and visualizing full-order solutions.

We can extract a mass matrix from this solver, but not the differential operators due to the velocity and source terms.

To animate the solution instead of redoing the computation:

```matlab
load('solver.mat', 'solver');
load('solution.mat', 'u');
solver.Animate_Solution(u.NodalSolution);
```

## Reduced-order Modeling

The script [Driver_Opt_Opinf.m](./Driver_Opt_OpInf.m) uses the `Transient_Adv_Diff_2D` solver to generate training data, learn a [POD basis](../../src/model_reduction/POD_Basis.m), and calibrate a bilinear reduced-order model through [Operator Inference](../../src/model_reduction/OpInf_ROM_Constraint.m):

$$
\begin{aligned}
    \frac{\textrm{d}}{\textrm{d}t}\hat{\mathbf{y}}(t)
    = \hat{\mathbf{A}}\hat{\mathbf{y}}(t) + \hat{\mathbf{B}}\mathbf{q}(t),
    \qquad\qquad
    \mathbf{y}(t) \approx \mathbf{V}_{\!r}\hat{\mathbf{y}}(t).
\end{aligned}
$$

The regression is populated with sixth-order finite difference estimates of the time derivatives of the training states.
