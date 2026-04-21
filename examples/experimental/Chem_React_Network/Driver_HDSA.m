%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));
rng(1423);

plot_figures = false;

T = 30;
n_t = 100;
con = React_Rate_Eqn(T, n_t);
obj = Chem_React_Network_Objective(T, n_t, con);
opt = Reduced_Space_Optimization(obj, con);

%%
data_interface = MD_Data_Interface_Chem_React_Network();
data_interface.Load_Data();

alpha_u = 0.0042;
alpha_z = 0.0027;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Chem_React_Network(alpha_u, opt);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Chem_React_Network(alpha_z);

%%
num_prior_samples = 500;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

prior_delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
t = linspace(0, T, n_t)';
if plot_figures
    figure;
    hold on;
    plot(t, prior_delta_samples(1:9:end, :), 'LineWidth', 3, 'color', [.9, .9, .9]);
    plot(t, prior_delta_samples(1:9:end, 1:5), 'LineWidth', 3);
    xlabel('Time');
    ylabel('$x_1$', 'Interpreter', 'latex');
    set(gca, 'fontsize', 18);
end

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
num_post_samples = 500;
alpha_d = 5.2405e-10;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(1, 2);
Z_test(:, 1) = data_interface.z_opt;
Z_test(:, 2) = 3;
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

if plot_figures
    for j = 1:9
        figure;
        hold on;
        plot(t, md_post_sampling.post_data.D(j:9:end, 1), 'color', 'black', 'LineWidth', 3);
        plot(t, delta_mean{1}(j:9:end), '--', 'color', 'red', 'LineWidth', 3);
        for k = 1:num_post_samples
            plot(t, delta_samples{1}(j:9:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
        end
        plot(t, md_post_sampling.post_data.D(j:9:end, 1), 'color', 'black', 'LineWidth', 3);
        plot(t, delta_mean{1}(j:9:end), '--', 'color', 'red', 'LineWidth', 3);
        xlabel('Time');
        ylabel(['Species ', num2str(j)], 'Interpreter', 'latex');
        set(gca, 'fontsize', 18);
    end
end

%%
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);
[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

%%
SSA = SSA_System(con);
num_samples = 500;
u_SSA_opt = SSA.SSA_Mean(data_interface.z_opt, num_samples);
y_SSA_opt = reshape(u_SSA_opt, 9, n_t)' * con.nA * con.vol / con.state_scale;
u_SSA_update_mean = SSA.SSA_Mean(z_update_mean, num_samples);
y_SSA_update_mean = reshape(u_SSA_update_mean, 9, n_t)' * con.nA * con.vol / con.state_scale;

disp('Error in terminal species with nominal:');
disp(abs(y_SSA_opt(end, 5) - (obj.target * con.nA * con.vol / con.state_scale)) / (obj.target * con.nA * con.vol / con.state_scale));

disp('Error in terminal species with update:');
disp(abs(y_SSA_update_mean(end, 5) - (obj.target * con.nA * con.vol / con.state_scale)) / (obj.target * con.nA * con.vol / con.state_scale));
