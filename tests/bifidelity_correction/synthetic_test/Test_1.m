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

z_update_ref = load('reference_solution.mat').z_update;
ref_diff = norm(z_update - z_update_ref);
if ref_diff > 1.e-9
    fprintf(2, '\nbifidelity_correction/synthetic_test failed.\n');
else
    fprintf(1, '\nbifidelity_correction/synthetic_test passed.\n');
end

% save('reference_solution.mat','z_update')

% figure;
% hold on;
% plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3);
% plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
% plot(x, z_update, '--', 'color', 'red', 'LineWidth', 3);