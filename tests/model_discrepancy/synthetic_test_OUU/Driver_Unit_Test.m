%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

suppress_figures = false;

m = 51;

data_interface = MD_OUU_Data_Interface_synthetic_test_OUU();
data_interface.Load_Data();

Xi = data_interface.Xi;
N = size(Xi, 2);
obj = Synthetic_Test_OUU_Objective(m);
cons = cell(N, 1);
for k = 1:N
    cons{k} = Synthetic_Test_OUU_Constraint(Xi(:, k));
end
opt = Reduced_Space_Optimization_Under_Uncertainty(obj, cons);
x = obj.x;

opt_prob_interface = MD_OUU_Opt_Prob_Interface_Sabl(data_interface, opt);

us_prior_interface = MD_u_Prior_Interface_synthetic_test_OUU(m);
u_prior_interface = MD_OUU_u_Prior_Interface(us_prior_interface, data_interface);
z_prior_interface = MD_z_Prior_Interface_synthetic_test_OUU(m);

%%
error = [];
M_u = kron(us_prior_interface.M,u_prior_interface.C);
W_u = kron(us_prior_interface.W_u,u_prior_interface.C);

%%
u_test = randn(m*N,1);
local_error = norm(M_u * u_test - u_prior_interface.Apply_M_u(u_test))/norm(M_u * u_test)
error = [error ; local_error];

%%
local_error = norm(linsolve(W_u,u_test) - u_prior_interface.Apply_W_u_Inverse(u_test))/norm(linsolve(W_u,u_test))
error = [error ; local_error];

%%
scalar = randn;
local_error = norm(linsolve(W_u+scalar*M_u,u_test) - u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u_test,scalar))/norm(linsolve(W_u+scalar*M_u,u_test))
error = [error ; local_error];

%%
S = 10000;
u_samples = u_prior_interface.Sample_with_Covariance_W_u_Inverse(S);
cov_approx = cov(u_samples');
local_error = norm(W_u * cov_approx - eye(m*N))/norm(W_u)
error = [error ; local_error];

%%
u_samples = u_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(S,scalar);
cov_approx = cov(u_samples');
local_error = norm((W_u + scalar*M_u) * cov_approx - eye(m*N))/norm((W_u + scalar*M_u))
error = [error ; local_error];

%%
z = data_interface.z_opt;
[~,~,hessian_data] = opt.Jhat(z);
v = randn(m,1);
Hv = opt.Jhat_hessVec(hessian_data,v);
local_error = norm(opt_prob_interface.Apply_RS_Hessian(v,z) - Hv)/norm(Hv)
error = [error ; local_error];

