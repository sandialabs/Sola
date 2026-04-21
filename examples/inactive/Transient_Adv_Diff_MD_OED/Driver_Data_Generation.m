%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));

datafile = 'OpInf_Training_Data.mat';
show_states = false;

%% Set up the optimization problem.
n_y = 200;
n_t = 51;
T = 1;
num_space_control_nodes = 10;
n_z = num_space_control_nodes * (n_t - 1);
con = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes);

%% Solve the state equation for several random controls.
num_solves = 3;
Z = randn(n_z, num_solves);
Y = zeros(n_y * n_t, num_solves);
for k = 1:num_solves
    Y(:, k) = con.State_Solve(Z(:, k));
end

%% Visualize each set of drawn states.
if show_states
    x = con.x;
    t = con.t_mesh;
    for j = 1:num_solves
        u_reshape = reshape(Y(:, j), n_y, n_t);
        figure;
        for k = 1:n_t
            plot(x, u_reshape(:, k), 'LineWidth', 3);
            ylim([-.1 .1]);
            title(['Forward solve ', num2str(j), ' at time ', num2str(t(k))]);
            pause(.05);
        end
    end
end

%% Save the data.
save(datafile, "n_y", "n_z", "T", "n_t", "num_space_control_nodes", "num_solves", "Y", "Z");
