clear;
close all;
clc;
addpath(genpath('../../../src/'));
% rng(154);

suppress_figures = false;

load Optimization_Results.mat;

m = 100;
theta = 1 + (1 - linspace(0, 1, m)').^2;
con = Adv_Diff_Constraint(theta);
% Generate_Obs_Data(con);
likelihood = Adv_Diff_Likelihood_Model(m);
prior = Adv_Diff_Prior_Model(con);
bayes_inv = Bayesian_Inversion(likelihood, prior, con);

gevp_fd = Randomized_GEVP_FD_Hessian(bayes_inv.opt, z_opt);
gevp_full = Randomized_GEVP_Full_Hessian(bayes_inv.opt, z_opt);

num_evals = 20;
oversampling = 20;
[evecs_fd, evals_fd] = gevp_fd.Compute_GEVP(num_evals, oversampling);
[evecs_full, evals_full] = gevp_full.Compute_GEVP(num_evals, oversampling);

H = gevp_full.Apply_Operator(eye(m));

[V, D] = eig(H, gevp_fd.M, 'vector');
[~, ind] = sort(D, 'descend');
D = D(ind);
V = V(:, ind);
V = V(:, 1:num_evals);
D = D(1:num_evals);
for k = 1:num_evals
    V(:, k) = V(:, k) / sqrt(V(:, k)' * gevp_fd.M * V(:, k));
end

low_rank_approx_error_full = norm(gevp_fd.M * V * diag(D) * V' * gevp_fd.M - gevp_fd.M * evecs_full * diag(evals_full) * evecs_full' * gevp_fd.M, 'fro') / norm(gevp_fd.M * V * diag(D) * V' * gevp_fd.M, 'fro');
low_rank_approx_error_fd = norm(gevp_fd.M * V * diag(D) * V' * gevp_fd.M - gevp_fd.M * evecs_fd * diag(evals_fd) * evecs_fd' * gevp_fd.M, 'fro') / norm(gevp_fd.M * V * diag(D) * V' * gevp_fd.M, 'fro');

disp('Low rank approximation error in exact hess-vec:');
disp(low_rank_approx_error_full);
disp('Low rank approximation error in finite difference hess-vec:');
disp(low_rank_approx_error_fd);

% H_FD = gevp_fd.Apply_Operator(eye(m));
% hessian_fd_error = norm(H-H_FD,'fro')/norm(H,'fro');
% H_FD_Sym = 0.5 * (H_FD + H_FD');
% hessian_fd_sym_error = norm(H-H_FD_Sym,'fro')/norm(H,'fro');

% [V_FD_Sym, D_FD_Sym] = eig(H_FD_Sym, gevp_fd.M, 'vector');
% [~,ind] = sort(D_FD_Sym,'descend');
% D_FD_Sym = D_FD_Sym(ind);
% V_FD_Sym = V_FD_Sym(:,ind);
% V_FD_Sym = V_FD_Sym(:, 1:num_evals);
% D_FD_Sym = D_FD_Sym(1:num_evals);
% for k = 1:num_evals
%     V_FD_Sym(:,k) = V_FD_Sym(:,k)/sqrt(V_FD_Sym(:,k)'*gevp_fd.M*V_FD_Sym(:,k));
% end
% low_rank_approx_error_FD_Sym = norm(gevp_fd.M*V*diag(D)*V'*gevp_fd.M - gevp_fd.M*V_FD_Sym*diag(D_FD_Sym)*V_FD_Sym'*gevp_fd.M,'fro')/norm(gevp_fd.M*V*diag(D)*V'*gevp_fd.M,'fro');
