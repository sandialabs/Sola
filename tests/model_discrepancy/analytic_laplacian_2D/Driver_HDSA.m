%% Set up
clear;
addpath(genpath('../../../src'));
load Optimization_Results.mat;
rng(2451423);

x = adv_diff.pde_meshing.x;
y = adv_diff.pde_meshing.y;
M = adv_diff.pde_meshing.M;
S = adv_diff.pde_meshing.S;
m = length(x);

data_interface = MD_Data_Interface_Adv_Diff();
data_interface.Load_Data();

hyperparams = MD_Hyperparameters_hyperparam_auto_2D(data_interface,x,y);
u_prior_interface = MD_Analytic_Laplacian_u_Prior_Interface(M,hyperparams);

hyperparams_gsvd = MD_Hyperparameters_hyperparam_auto_2D(data_interface,x,y);
hyperparams_gsvd.gsvd_num_sing_vals = 1849;
u_prior_interface_gsvd = MD_Laplacian_u_Prior_Interface(S,M,hyperparams_gsvd);


%%
v = randn(m,1);
u1 = u_prior_interface.Apply_W_u_Inverse(v);
u2 = u_prior_interface_gsvd.Apply_W_u_Inverse(v);
diff1 = norm(u1-u2)/norm(u1);

if diff1 > 1.e-2
    disp('analytic_laplacian_2D difference:');
    disp(diff1);
end

v = randn(m,1);
scalar = rand;
u1 = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(v,scalar);
u2 = u_prior_interface_gsvd.Apply_W_u_Plus_scalar_M_u_Inverse(v,scalar);
diff2 = norm(u1-u2)/norm(u1);

if diff2 > 1.e-2
    disp('analytic_laplacian_2D difference:');
    disp(diff2);
end

%%
i = 3;
j = 2;
u = cos(i*(pi/2)*(x+1)).*cos(j*(pi/2)*(y+1));
z = u*((1+u_prior_interface.beta_u*(pi/2)^2*(i^2+j^2))^2)/u_prior_interface.alpha_u;

v = u_prior_interface.Apply_M_u(z);
u_approx = u_prior_interface.Apply_W_u_Inverse(v);

diff3 = norm(u-u_approx,'fro')/norm(u,'fro');

if diff3 > 1.e-3
    disp('analytic_laplacian_2D difference:');
    disp(diff3);
end