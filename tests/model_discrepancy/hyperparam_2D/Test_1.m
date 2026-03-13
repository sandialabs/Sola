%% Set up
clear;
close all;

addpath(genpath('../../../src'));
load Optimization_Results.mat;
rng(2451423);

x = adv_diff.pde_meshing.x;
y = adv_diff.pde_meshing.y;
M = adv_diff.pde_meshing.M;
S = adv_diff.pde_meshing.S;
m = length(x);

obj = Adv_Diff_Objective(adv_diff, reg_coeff);
con = Adv_Diff_Constraint(adv_diff);
opt = Reduced_Space_Optimization(obj, con);

data_interface = MD_Data_Interface_Adv_Diff();
data_interface.Load_Data();

u_hyperparams = MD_u_Hyperparameter_Interface_hyperparam_2D(x, y);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparams);

% save('reference_solution.mat','u_hyperparams')
ref = load('reference_solution.mat');
ref_diff = norm(u_hyperparams.alpha_u - ref.u_hyperparams.alpha_u);
ref_diff = max(ref_diff, norm(u_hyperparams.beta_u - ref.u_hyperparams.beta_u));

if ref_diff > 1.e-9
    fprintf(2,'\nmodel_discrepancy/hyperparam_2D failed.\n');
else
    fprintf(1,'\nmodel_discrepancy/hyperparam_2D passed.\n');
end
