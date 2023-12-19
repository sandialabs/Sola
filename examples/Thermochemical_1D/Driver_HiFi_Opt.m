clear;
close all;
clc;
run('../../src/Set_Paths');

n_y = 50;
T = 1;
n_t = 51;
control_time_nodes = 10;
n_z = n_y * control_time_nodes;
reg_coeff = 1.e-3;

con = Thermochemical_HiFi_Constraint_AD(n_y, control_time_nodes, T, n_t);
con.AD_Initialization('Hifi_AD_Files');
obj = Thermochemical_Dynamic_Objective(con, reg_coeff);
opt = Reduced_Space_Optimization(obj, con);

z0 = 10 * ones(n_z, 1);
[u, z] = opt.Optimize(z0);

f = con.Map_Controller_to_Mesh(z);
u_rs = reshape(u, 6 * n_y, n_t);
T = con.I_T * u_rs;
u1 = con.I_u1 * u_rs;
u2 = con.I_u2 * u_rs;
v1 = con.I_v1 * u_rs;
v2 = con.I_v2 * u_rs;
v3 = con.I_v3 * u_rs;

figure(1);
surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', T);
xlabel('Space');
ylabel('Time');
title('T');
colorbar();
set(gca, 'fontsize', 18);

figure(2);
surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', u1);
xlabel('Space');
ylabel('Time');
title('u1');
colorbar();
set(gca, 'fontsize', 18);

figure(3);
surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', u2);
xlabel('Space');
ylabel('Time');
title('u2');
colorbar();
set(gca, 'fontsize', 18);

figure(4);
surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', v1);
xlabel('Space');
ylabel('Time');
title('v1');
colorbar();
set(gca, 'fontsize', 18);

figure(5);
surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', v2);
xlabel('Space');
ylabel('Time');
title('v2');
colorbar();
set(gca, 'fontsize', 18);

figure(6);
surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', v3);
xlabel('Space');
ylabel('Time');
title('v3');
colorbar();
set(gca, 'fontsize', 18);

figure(7);
surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', f);
xlabel('Space');
ylabel('Time');
title('f');
colorbar();
set(gca, 'fontsize', 18);

rmdir('Hifi_AD_Files/', 's');
save('HiFi_Opt_Results.mat', 'obj', 'con', 'u', 'z');
