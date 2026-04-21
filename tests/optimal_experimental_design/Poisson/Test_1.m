%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;
rng(1423435);

print_output = false;

m = 200;
diff_coeff = 1;
con = Poisson_Constraint(m, diff_coeff);

sigma = 1.e-1;
obs_vec = 9:33;
likelihood = Poisson_Likelihood_Model(sigma, obs_vec, m);
prior = Poisson_Prior_Model(con);

num_trace_samples = 1000;
reguarlization_coeff = 1.e-1;
linear_oed = Linear_OED(likelihood, prior, con, num_trace_samples, reguarlization_coeff);

max_error = 0;

%%
num_sing_vals = length(obs_vec);
oversampling = 0;
num_subspace_iters = 1;
sing_vals = linear_oed.Compute_Forward_Operator_GSVD(num_sing_vals, oversampling, num_subspace_iters);

U = linear_oed.forward_operator_sing_vecs_output;
Sigma = diag(linear_oed.forward_operator_sing_vals);
V = linear_oed.forward_operator_sing_vecs_input;
prior_precision = prior.Prior_Precision_Apply(eye(m));
forward_approx = U * Sigma * V' * prior_precision;

tmp = linsolve(con.A, con.B);
forward = tmp(obs_vec, :);
error = norm(forward - forward_approx);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
w = round(rand(length(obs_vec), 1));
[evecs, evals] = linear_oed.Compute_Misfit_Hessian_GEVP(w);
Minv = prior.Mass_Matrix_Inverse_Apply(eye(m));
rank = min(length(find(w ~= 0)), length(obs_vec));
error = norm(evecs' * Minv * evecs - eye(rank));
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
misfit_hessian_approx = Minv * evecs * diag(evals) * evecs' * Minv;

Ftilde = (1 / likelihood.sigma) * forward * prior.Laplacian_Like_Inverse_Apply(eye(m));
misfit_hessian = Ftilde' * diag(w) * Ftilde;
error = norm(misfit_hessian - misfit_hessian_approx);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
L = prior.Laplacian_Like_Apply(eye(m));
H_approx = L' * (Minv + misfit_hessian) * L;

likelihood.obs_vec = likelihood.obs_vec(find(w ~= 0));
bayes_inv = Bayesian_Inversion(likelihood, prior, con);
z0 = randn(m, 1);
[~, ~, hessian_data] = bayes_inv.opt.Jhat(z0);
H = bayes_inv.opt.Jhat_hessVec(hessian_data, eye(m));
error = norm(H - H_approx);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
v = randn(m, 1);
Hinv_v_approx = linear_oed.Compute_Inverse_Hessian_Matvec(v, evecs, evals);

Hinv_v = linsolve(H, v);
error = norm(Hinv_v - Hinv_v_approx);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
if print_output
    disp(['error = ', num2str(error)]);
end

%%
w = rand(length(obs_vec), 1);
[val, grad] = linear_oed.Posterior_Trace_Objective(w);
dw = randn(length(grad), 1);
dw = dw / (max(abs(dw)) + .5);
h = 10.^(-1:-1:-6)';
M = length(h);
fd = zeros(M, 1);
for k = 1:M
    val_k = linear_oed.Posterior_Trace_Objective(w + h(k) * dw);
    fd(k) = (val_k - val) / h(k);
end
if print_output
    disp('Finite difference step sizes (log10)');
    disp(log10(h)');
    disp('Finite difference errors (log10)');
    disp(log10(abs(fd - grad' * dw)'));
end

%%
w = rand(length(obs_vec), 1);
[val, grad] = linear_oed.OED_Objective(w);
dw = randn(length(grad), 1);
dw = dw / (max(abs(dw)) + .5);
h = 10.^(-1:-1:-6)';
M = length(h);
fd = zeros(M, 1);
for k = 1:M
    val_k = linear_oed.OED_Objective(w + h(k) * dw);
    fd(k) = (val_k - val) / h(k);
end
if print_output
    disp('Finite difference step sizes');
    disp(log10(h)');
    disp('Finite difference errors (log10)');
    disp(log10(abs(fd - grad' * dw)'));
end

if error > 1.e-9
    fprintf(2, '\noptimal_experimental_design/Poisson failed.\n');
else
    fprintf(1, '\noptimal_experimental_design/Poisson passed.\n');
end
