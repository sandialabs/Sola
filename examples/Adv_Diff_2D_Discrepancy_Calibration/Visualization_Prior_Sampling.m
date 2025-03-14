%% Set up
function [x,y, delta_samples_z_opt, delta_samples_z_pert, z_pert, corlength, mags] = Visualization_Prior_Sampling()
addpath(genpath('../../src'));
load Optimization_Results.mat;
load Assembled_Operators.mat;
x = pde_meshing.x;
y = pde_meshing.y;
m = length(x);

con = Diff_Constraint(pde_meshing, diff_coeff);
obj = Diff_Objective(con, reg_coeff);

data_interface = MD_Data_Interface_Diff();
data_centering = true;

u_hyperparams = MD_u_Hyperparameters_Diff(data_interface, x, y, data_centering);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(pde_meshing.S, pde_meshing.M, u_hyperparams);

num_state_solves = 100;
z_hyperparams = MD_z_Hyperparameters_Diff(data_interface, u_prior_interface, num_state_solves, x, y, con);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(pde_meshing.S, pde_meshing.M, z_hyperparams);

z_hyperparams.beta_z = (5) * z_hyperparams.beta_z;
z_prior_interface.Set_beta_z(z_hyperparams.beta_z);

z_hyperparams.alpha_z = (1/4) * z_hyperparams.alpha_z;
z_prior_interface.Set_alpha_z(z_hyperparams.alpha_z);

md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

num_prior_samples = 100;
num_perts_init = round((18/pi^2)/z_prior_interface.beta_z);
[delta_samples_z_opt, delta_samples_z_pert, z_pert] = md_prior_sampling.Prior_Discrepancy_Samples_for_Visualization(num_prior_samples, num_perts_init);
num_perts = size(z_pert,2);

corlength = zeros(num_perts,1);
mags = zeros(num_prior_samples, num_perts);

initial_guess = 0;
for i = 1:num_perts
    corlength(i) = computeCorrelationLength_2D(x,y,z_pert(:,i), initial_guess);
    initial_guess = corlength(i);
    fprintf("%d : %.3f\n", i, corlength(i));

    for j = 1:num_prior_samples
        mags(j,i) = max(abs(delta_samples_z_pert{i}(:,j)));
    end
end

corlength_raw = corlength;
L = 18;
e = linspace(min(corlength)-max(corlength)/40,max(corlength)+max(corlength)/40,L)';
for i = 1:num_perts
    I = find(e < corlength(i));
    k = I(end);
    corlength(i) = .5*(e(k)+e(k+1));
end

U = length(unique(corlength));
while (U < 18) && ((e(2)-e(1)) > (max(corlength_raw)-min(corlength_raw))/30)
    L = L + 1;
    corlength = corlength_raw;
    e = linspace(min(corlength)-max(corlength)/40,max(corlength)+max(corlength)/40,L)';
    for i = 1:num_perts
        I = find(e < corlength(i));
        k = I(end);
        corlength(i) = .5*(e(k)+e(k+1));
    end
    U = length(unique(corlength));
end

% name = 'Discrepancy sample 1 at z_{opt}';
% pde_meshing.Plot_Field(delta_samples_z_opt(:, 1), name);
% 
% name = 'Discrepancy sample 1 at pertubed z';
% pde_meshing.Plot_Field(delta_samples_z_pert{1}(:, 1), name);
% 
% name = 'Perturbed z';
% pde_meshing.Plot_Field(z_pert(:, 1), name);
