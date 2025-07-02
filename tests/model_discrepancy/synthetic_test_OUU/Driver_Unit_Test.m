%%
clear;
close all;
clc;
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
scalar = rand;
local_error = norm(linsolve(W_u+scalar*M_u,u_test) - u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u_test,scalar))/norm(linsolve(W_u+scalar*M_u,u_test))
error = [error ; local_error];

%%
theta_test = randn(m*(m+1)*N,1);
W_theta = zeros(m*(m+1),m*(m+1));
W_theta(1:m,1:m) = us_prior_interface.W_u;
v = us_prior_interface.M*data_interface.z_opt;
W_theta(1:m,(m+1):end) = kron(us_prior_interface.W_u,v');
W_theta((m+1):end,1:m) = kron(us_prior_interface.W_u,v);
W_theta((m+1):end,(m+1):end) = kron(us_prior_interface.W_u,z_prior_interface.W_z + v*v');

test1 = kron(W_theta,u_prior_interface.C) * theta_test;

theta_test_rs = reshape(theta_test,N,m*(m+1))';
test2 = u_prior_interface.C * theta_test_rs' * W_theta;
test2 = test2(:);
local_error = norm(test1 - test2)/norm(test1)
error = [error ; local_error];

E = zeros(N,N);
for i = 1:N
    for j = 1:N
        if i==j
            E(i,j) = (1/u_prior_interface.scaling) * u_prior_interface.L(i,j) * theta_test_rs(:,i)'*W_theta*theta_test_rs(:,j);
        else
            E(i,j) = (1/u_prior_interface.scaling) * u_prior_interface.L(i,j) * (theta_test_rs(:,i) - theta_test_rs(:,j))' * W_theta * (theta_test_rs(:,i) - theta_test_rs(:,j));
        end
    end
end
test3 = sum(E(:));
local_error = abs(theta_test'*test1 - test3)/abs(theta_test'*test1)
error = [error ; local_error];

%%
z = data_interface.z_opt;
[~,~,hessian_data] = opt.Jhat(z);
v = randn(m,1);
Hv = opt.Jhat_hessVec(hessian_data,v);
local_error = norm(opt_prob_interface.Apply_RS_Hessian(v,z) - Hv)/norm(Hv)
error = [error ; local_error];

%%
S = zeros(m*N,m);
for k = 1:N
    I = ((k-1)*m + 1):(k*m);
    S(I,:) = diag(Xi(1,k) * 3 * z.^2);
end
u_test = randn(m*N,1);
tmp = data_interface.Reshape_State_to_Mat(u_test);
test1 = S'*tmp(:);
test2 = opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_test,z);
local_error = norm(test1 - test2)/norm(test1)
error = [error ; local_error];

%%
u = data_interface.u_opt;
u_test = randn(m*N,1);
tmp = data_interface.Reshape_State_to_Mat(u_test);
u_out = zeros(m,N);
for k = 1:N
    u_out(:,k) = (1/N) * us_prior_interface.M * tmp(:,k);
end
test1 = data_interface.Reshape_State_to_Vec(u_out);
test2 = opt_prob_interface.Apply_Misfit_Hessian(u_test,u,z);
local_error = norm(test1 - test2)/norm(test1)
error = [error ; local_error];

%%
u = data_interface.u_opt;
u_tmp = data_interface.Reshape_State_to_Mat(u);
u_out = zeros(m,N);
for k = 1:N
    u_out(:,k) = (1/N) * us_prior_interface.M * (u_tmp(:,k) - obj.T);
end
test1 = data_interface.Reshape_State_to_Vec(u_out);
test2 = opt_prob_interface.Misfit_Gradient(u,z);
local_error = norm(test1 - test2)/norm(test1)
error = [error ; local_error];

disp('The following tests are sample covariance estimators. We should only expect local errors on the order of 1.e-2 or 1.e-3')
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