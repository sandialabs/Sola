# Transient Advection-Diffusion in Two Spatial Dimensions

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

Show geometry and edge labels: `pdegplot(model,EdgeLabels="on")`.

To animate the solution instead of redoing the computation:

```matlab
load('solver.mat', 'solver');
load('solution.mat', 'u');
solver.Animate_Solution(u.NodalSolution);
```
