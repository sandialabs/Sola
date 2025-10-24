%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

random_numbers = randn(3 * 10^6, 1);
writematrix(random_numbers, 'random_numbers.txt');

rng(121234);

m = 51;

data_interface = MD_OUU_Data_Interface_synthetic_test_OUU();
data_interface.Load_Data();

fileID = fopen('z_opt.txt', 'w');
fprintf(fileID, '%f %f\n', data_interface.z_opt');
fclose(fileID);

fileID = fopen('Z.txt', 'w');
fprintf(fileID, '%f %f\n', data_interface.Z');
fclose(fileID);

Xi = load('Optimization_Results.mat','Xi').Xi;
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
ensemble_weighting = MD_OUU_Ensemble_Weighting_Matrix(data_interface, us_prior_interface);
ensemble_weighting.Set_Matrices(0.357);
u_prior_interface = MD_OUU_u_Prior_Interface(us_prior_interface, data_interface, ensemble_weighting);
z_prior_interface = MD_z_Prior_Interface_synthetic_test_OUU(m);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
prior_delta_z_opt = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

%%
z = zeros(m, 3);
z(:, 1) = x;
z(:, 2) = x.^2 + 1;
z(:, 3) = sin(2 * pi * x);
prior_delta = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = randn(m, 3);
Z_test(:, 1:2) = md_post_sampling.post_data.Z;
Z_test(:, 3) = 1.5 * ones(m, 1);
[post_delta_mean, post_delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 8;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

md_update = MD_Update(md_post_sampling, md_hessian_analysis);
[post_z_mean, post_z_samples] = md_update.Posterior_Update_Samples();

%%
post_delta_mean = reshape(cell2mat(post_delta_mean), 1530, 3);
save('Sabl_Output.mat', 'prior_delta', 'prior_delta_z_opt', 'post_delta_mean', 'post_delta_samples', 'post_z_mean', 'post_z_samples');
