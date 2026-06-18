%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;
rng(121234);

m = 51;
x = linspace(0, 1, m)';
z_opt = 1 + x;

opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test(m);
sol_op_interface = BF_Sol_Op_Interface_synthetic_test();

bf_update = BF_Update(sol_op_interface, opt_prob_interface);

z_update = bf_update.Update(z_opt);

figure;
hold on;
plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3);
plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update, '--', 'color', 'red', 'LineWidth', 3);