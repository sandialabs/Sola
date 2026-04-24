%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;
rng(121234);

suppress_figures = true;
error = [];

m = 51;
x = linspace(0, 1, m)';

[M, S] = Assemble_Mass_and_Stiffness(m);
M = sparse(M);
S = sparse(S);

data_interface = MD_Data_Interface_synthetic_test_bilaplacian(m);

u_hyperparam_interface = MD_u_Hyperparameter_Interface_synthetic_test_bilaplacian(m);
u_hyperparam_interface.alpha_u = 0.048969233204560;
u_hyperparam_interface.beta_u = 0.007702351792463;
u_hyperparam_interface.alpha_d = 2.177109166165424e-07;
u_prior_interface = MD_Bilaplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

scalar = 0.5;
u_in = linspace(0, 1, m)';
num_samples = 1;

E_u = u_prior_interface.E_u;
W_u_acute = E_u' * (M \ eye(m)) * E_u;

u_out = u_prior_interface.Apply_M_u(u_in);
test = M * u_in;
local_error = norm(u_out - test) / norm(test);
error = [error; local_error];

u_out = u_prior_interface.Apply_W_u_Acute_Plus_scalar_M_u_Inverse(u_in, scalar);
test = (W_u_acute + scalar * M) \ u_in;
local_error = norm(u_out - test) / norm(test);
error = [error; local_error];

u_out = u_prior_interface.Apply_W_u_Acute_Inverse(u_in);
test = W_u_acute \ u_in;
local_error = norm(u_out - test) / norm(test);
error = [error; local_error];

u_sample_1 = u_prior_interface.Sample_with_Covariance_W_u_Acute_Inverse(num_samples);

u_sample_2 = u_prior_interface.Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(num_samples, scalar);

rng(121234);

omega = randn(m, num_samples);
L = u_prior_interface.R_mass';
Linv = L \ eye(m);
S = L * sqrtm(full(Linv * M * Linv'));
vec = S * omega;
u_sample_1_test = E_u \ vec;
local_error = norm(u_sample_1 - u_sample_1_test) / norm(u_sample_1_test);
error = [error; local_error];

omega = randn(m, num_samples);
A = W_u_acute + scalar * M;
A_approx = u_prior_interface.W_u_acute_approx + scalar * M;
L = ichol(A_approx);
Linv = L \ eye(m);
S = L * sqrtm(full(Linv * A * Linv'));
vec = S * omega;
u_sample_2_test = A \ vec;

local_error = norm(u_sample_2 - u_sample_2_test) / norm(u_sample_2_test);
error = [error; local_error];

if max(error) > 1.e-7
    fprintf(2, '\nmodel_discrepancy/synthetic_test_bilaplacian failed.\n');
else
    fprintf(1, '\nmodel_discrepancy/synthetic_test_bilaplacian passed.\n');
end
