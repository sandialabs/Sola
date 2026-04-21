%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;

n_y = 31;
n_t = 12;

T = 1;
x = linspace(0, 1, n_y)';
t = linspace(0, T, n_t)';

[M_s, S_s] = Assemble_Mass_and_Stiffness(n_y);
[M_t, S_t] = Assemble_Mass_and_Stiffness(n_t);

data_interface = MD_Data_Interface_lumped_mass_unit_test(n_y, n_t);

u_hyperparam_interface = MD_u_Hyperparameter_Interface_synthetic_test_lumped_mass(n_y);
u_hyperparam_interface.alpha_u = 0.048969233204560;
u_hyperparam_interface.beta_u = 0.007702351792463;
u_hyperparam_interface.beta_t = 0.1;

spatial_u_prior_interface = MD_Lumped_Mass_u_Prior_Interface(S_s, M_s, data_interface, u_hyperparam_interface);
transient_prior_cov = MD_Transient_Prior_Covariance_Sola(data_interface, u_hyperparam_interface, T, n_t, n_y);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(data_interface, spatial_u_prior_interface, transient_prior_cov);

%%
E_s = full(u_hyperparam_interface.beta_u * S_s + M_s);
E_s_inv = linsolve(E_s, eye(n_y));

M_lumped = diag(M_s * ones(n_y, 1));
M_lumped_inv = linsolve(M_lumped, eye(n_y));
W_s = (1 / u_hyperparam_interface.alpha_u) * E_s' * M_lumped_inv * E_s;
W_s_inv = u_hyperparam_interface.alpha_u * E_s_inv * M_lumped * E_s_inv';

W_t = full(u_hyperparam_interface.beta_t * S_t + M_t);
W_t_inv = linsolve(W_t, eye(n_t));
M_t = full(M_t);

M_u = kron(M_t, M_s);
W_u = kron(W_t, W_s);
W_u_inv = kron(W_t_inv, W_s_inv);

diff = [];

%%
Lambda = diag(1 ./ transient_prior_cov.evals);
V = transient_prior_cov.evecs;

test = M_t * V * Lambda * V' * M_t;
local_diff = norm(test - W_t, 'fro') / norm(W_t, 'fro');
diff = [diff; local_diff];

%% Apply_M_u
u = x .* t';
u = u(:);
test = u' * kron(M_t, M_s) * u;
local_diff = abs(test - 1 / 9);
diff = [diff; local_diff];

test = u' * u_prior_interface.Apply_M_u(u);
local_diff = abs(test - 1 / 9);
diff = [diff; local_diff];

%% Apply_W_u_Inverse
u = randn(n_y * n_t, 1);
test1 = u_prior_interface.Apply_W_u_Inverse(u);
test2 = W_u_inv * u;
local_diff = norm(test1 - test2) / norm(test2);
diff = [diff; local_diff];

%% Apply_W_u_Plus_scalar_M_u_Inverse
u = randn(n_y * n_t, 1);
beta = (1.e5) * rand;
test1 = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u, beta);
test2 = linsolve(W_u + beta * M_u, u);
local_diff = norm(test1 - test2) / norm(test2);
diff = [diff; local_diff];

%% Sample_with_Covariance_W_u_Inverse

%%
T = kron(V * inv(sqrt(Lambda)), eye(n_y));
local_diff = norm(T * T' - kron(W_t_inv, eye(n_y)), 'fro') / norm(kron(W_t_inv, eye(n_y)), 'fro');
diff = [diff; local_diff];

%%
L = sqrt(u_hyperparam_interface.alpha_u) * kron(eye(n_t), E_s_inv * sqrt(M_lumped));
S = T * L;

local_diff = norm(S * S' - W_u_inv, 'fro') / norm(W_u_inv, 'fro');
diff = [diff; local_diff];

%%
seed = randi(1.e5);
rng(seed);
num_samples = 2;
u_samples = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samples);

rng(seed);
omega1 = randn(n_t * n_y, 1);
omega2 = randn(n_t * n_y, 1);

local_diff1 = norm(S * omega1 - u_samples(:, 1), 'fro') / norm(S * omega1, 'fro');
local_diff2 = norm(S * omega2 - u_samples(:, 2), 'fro') / norm(S * omega2, 'fro');
diff = [diff; local_diff1];
diff = [diff; local_diff2];

%% Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse

%%
beta = (1.e5) * rand;
A = sparse(W_s + beta * M_s);
L = ichol(A);
W_u_Acute_Plus_scalar_M_u_sqrt = Sparse_Matrix_Sqrt(A, L);
omega = randn(n_y, 1);
u_test = W_u_Acute_Plus_scalar_M_u_sqrt.Matrix_Sqrt_Apply(omega);

G = L \ eye(n_y);
B = G * A * G';
u = G \ sqrtm(B) * omega;

local_diff = norm(u - u_test) / norm(u);
diff = [diff; local_diff];

%%
L = sqrtm(linsolve(kron(eye(n_t), W_s) + kron(beta * inv(Lambda), M_s), eye(n_t * n_y)));
S = T * L;
W_u_beta_M_u_inv = linsolve(W_u + beta * M_u, eye(n_t * n_y));

local_diff = norm(S * S' - W_u_beta_M_u_inv, 'fro') / norm(W_u_beta_M_u_inv, 'fro');
diff = [diff; local_diff];

%%
seed = randi(1.e5);
rng(seed);
num_samples = 1;
u_prior_interface.spatial_prior_cov.use_sampling_prec = false;
u_samples = u_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(num_samples, beta);

rng(seed);
omega = randn(n_t * n_y, 1);

local_diff = norm(S * omega - u_samples(:, 1), 'fro') / norm(S * omega, 'fro');
diff = [diff; local_diff];

%%
if max(diff) > 1.e-8
    fprintf(2, '\nmodel_discrepancy/lumped_mass_unit_test failed.\n');
else
    fprintf(1, '\nmodel_discrepancy/lumped_mass_unit_test passed.\n');
end
