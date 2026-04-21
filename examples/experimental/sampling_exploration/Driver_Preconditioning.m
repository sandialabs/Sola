%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

M = load('fem_matrices.mat', 'mass_matrix').mass_matrix;
S = load('fem_matrices.mat', 'stiffness_matrix').stiffness_matrix;

% [M, S] = Assemble_Mass_and_Stiffness(100);

E = (0.01) * S + M;
m = size(M, 1);
M_lump = M * ones(m, 1);
beta = 1.e4;
W_lumped_mat = E * sparse(diag(1 ./ M_lump)) * E + beta * M;
tol = 1.e-8;

L = ichol(W_lumped_mat, struct('type', 'ict', 'droptol', 1e-7));

%%
b = randn(m, 1);

W = @(x) E * M \ (E * x) + beta * M * x;
W_lumped = @(x) W_lumped_mat * x;

[x_W, relres_W] = krylov_sqrt(W, b, m, tol);

A = @(x) L \ (W(L' \ x));
[x_W_pre, relres_W_pre] = krylov_sqrt(A, b, m, tol);
x_W_pre = L * x_W_pre;

Wlump_sqrt = Sparse_Matrix_Sqrt(W_lumped_mat);
[x_Wlump, relres_Wlump] = Wlump_sqrt.Matrix_Sqrt_Apply(b);

Wlump_pre_sqrt = Sparse_Matrix_Sqrt(W_lumped_mat, L);
[x_Wlump_pre, relres_Wlump_pre] = Wlump_pre_sqrt.Matrix_Sqrt_Apply(b);

disp('Number of iterations for:');
disp(['W = ', num2str(size(relres_W, 1))]);

disp('Number of iterations for:');
disp(['W_pre = ', num2str(size(relres_W_pre, 1))]);

disp('Number of iterations for:');
disp(['W_lump = ', num2str(size(relres_Wlump{1}, 1))]);

disp('Number of iterations for:');
disp(['Wlump_pre = ', num2str(size(relres_Wlump_pre{1}, 1))]);
