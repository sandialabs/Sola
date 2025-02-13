%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

suppress_figures = true;

m = 51;
x = linspace(0, 1, m)';

[M,S] = Assemble_Mass_and_Stiffness(m);

data_interface = MD_Data_Interface_synthetic_test_with_hyperparam_auto(m);
data_interface.Load_Data();

hyperparams = MD_Hyperparameters_synthetic_test_with_hyperparam_auto(data_interface,m);

u_prior_interface_gsvd = MD_Numeric_Laplacian_u_Prior_Interface(S,M,hyperparams);
u_prior_interface = MD_Analytic_Laplacian_u_Prior_Interface(M,hyperparams);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S,M,hyperparams);

%%
v = randn(m,1);
u1 = u_prior_interface.Apply_W_u_Inverse(v);
u2 = u_prior_interface_gsvd.Apply_W_u_Inverse(v);
diff1 = norm(u1-u2)/norm(u1);

if diff1 > 1.e-3
    disp('model_discrepancy_synthetic_test_with_analytic_laplacian difference:');
    disp(diff1);
end

v = randn(m,1);
scalar = rand;
u1 = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(v,scalar);
u2 = u_prior_interface_gsvd.Apply_W_u_Plus_scalar_M_u_Inverse(v,scalar);
diff2 = norm(u1-u2)/norm(u1);

if diff2 > 1.e-3
    disp('model_discrepancy_synthetic_test_with_analytic_laplacian difference:');
    disp(diff2);
end