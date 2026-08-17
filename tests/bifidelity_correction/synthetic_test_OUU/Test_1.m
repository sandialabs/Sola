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

data_interface = MD_OUU_Data_Interface_synthetic_test_OUU();
data_interface.Load_Data();
z_opt = data_interface.Load_Optimal_z();

Xi = load('Optimization_Results.mat', 'Xi').Xi;
M = size(Xi, 2);
obj = Synthetic_Test_OUU_Objective(m);
cons = cell(M, 1);
for k = 1:M
    cons{k} = Synthetic_Test_OUU_Constraint(Xi(:, k));
end
opt = Reduced_Space_Optimization_Under_Uncertainty(obj, cons);

ouu_opt_prob_interface = MD_OUU_Opt_Prob_Interface_Sola(data_interface, opt);
ouu_sol_op_interface = BF_OUU_Sol_Op_Interface_synthetic_test(data_interface.Load_Xi());

bf_update = BF_OUU_Update(ouu_sol_op_interface, ouu_opt_prob_interface);

z_update = bf_update.Update(z_opt);

z_update_ref = load('reference_solution.mat').z_update;
ref_diff = norm(z_update - z_update_ref);
if ref_diff > 1.e-9
    fprintf(2, '\nbifidelity_correction/synthetic_test_OUU failed.\n');
else
    fprintf(1, '\nbifidelity_correction/synthetic_test_OUU passed.\n');
end

% save('reference_solution.mat','z_update')

% z_hifi = load('Optimization_Results.mat','z_hifi').z_hifi;
% figure,
% plot(x,z_update,x,z_hifi,'--')

