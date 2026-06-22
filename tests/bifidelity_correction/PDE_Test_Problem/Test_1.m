%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;
rng(1234423);

m = 200;
diff_coeff = 1;
vel_coeff = 1 / 2;
robin_coeff = 2;
reg_coeff = 10;
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(obj, con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_hifi.x;

%%
data_interface = MD_Data_Interface_PDE_Test_Problem();
data_interface.Load_Data();
z_opt = data_interface.Load_Optimal_z();

sol_op_interface = BF_Sol_Op_Interface_Sola(con_hifi);
opt_prob_interface = MD_Opt_Prob_Interface_Sola(opt_lofi, data_interface);

bf_update = BF_Update(sol_op_interface, opt_prob_interface);

z_update = bf_update.Update(z_opt);

load z_hifi.mat
figure;
hold on;
plot(x, z_hifi, 'color', 'black', 'LineWidth', 3);
plot(x, z_opt, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update, '--', 'color', 'red', 'LineWidth', 3);

