clear;
clc;
close all;

M = load('fem_matrices.mat', 'mass_matrix').mass_matrix;
S = load('fem_matrices.mat', 'stiffness_matrix').stiffness_matrix;

% [M, S] = Assemble_Mass_and_Stiffness(50);

E = (0.01) * S + M;
m = size(M, 1);
M_lump = M * ones(m, 1);
beta = 1.e4;
W = E * sparse(diag(1 ./ M_lump)) * E + beta * M;
tol = 1.e-7;

%%
b = randn(m, 1);

A = @(x) W * x;
[x_W, relres_W] = krylov_sqrt(A, b, 1000, tol);

A = @(x) W \ x;
[x_Winv, relres_Winv] = krylov_sqrt(A, b, 1000, tol);

%%
M = 100;
omega = randn(m, M);

A = @(x) W * x;
x = zeros(m, M);
for k = 1:M
    x(:, k) = krylov_sqrt(A, omega(:, k), 1000, tol);
end
sample_1 = W \ x;

A = @(x) W \ x;
sample_2 = zeros(m, M);
for k = 1:M
    sample_2(:, k) = krylov_sqrt(A, omega(:, k), 1000, tol);
end

R = chol(W);
sample_3 = R \ omega;

Sigma_1 = cov(sample_1');
Sigma_2 = cov(sample_2');
Sigma_3 = cov(sample_3');
diff = norm(Sigma_1 - Sigma_2, 'fro') / norm(Sigma_2, 'fro');
diff1 = norm(Sigma_1 - Sigma_3, 'fro') / norm(Sigma_3, 'fro');
diff2 = norm(Sigma_2 - Sigma_3, 'fro') / norm(Sigma_3, 'fro');
