clear
close all

M = load('fem_matrices.mat','mass_matrix').mass_matrix;
S = load('fem_matrices.mat','stiffness_matrix').stiffness_matrix;

%[M, S] = Assemble_Mass_and_Stiffness(100);

E = (0.01) * S + M;
m = size(M,1);
M_lump = M * ones(m,1);
beta = 1.e4;
W_lumped_mat = E * sparse(diag(1./M_lump)) * E + beta * M;
tol = 1.e-7;

R = ichol(W_lumped_mat,struct('type','ict','droptol',1e-7));

%%
b = randn(m,1);

W = @(x) E * M \ (E * x) + beta * M * x;
W_lumped = @(x) W_lumped_mat * x;

[x_W, relres_W] = krylov_sqrt(W, b, 1000, tol);

A = @(x) R \ ( W(R' \ x) );
[x_W_pre, relres_W_pre] = krylov_sqrt(A, b, 1000, tol);

[x_Wlump, relres_Wlump] = krylov_sqrt(W_lumped, b, 1000, tol);

A = @(x) R \ ( W_lumped(R' \ x) );
[x_Wlump_pre, relres_Wlump_pre] = krylov_sqrt(A, b, 1000, tol);

disp('Number of iterations for:')
disp(['W = ',num2str(size(relres_W,1))])

disp('Number of iterations for:')
disp(['W_pre = ',num2str(size(relres_W_pre,1))])

disp('Number of iterations for:')
disp(['W_lump = ',num2str(size(relres_Wlump ,1))])

disp('Number of iterations for:')
disp(['Wlump_pre = ',num2str(size(relres_Wlump_pre ,1))])