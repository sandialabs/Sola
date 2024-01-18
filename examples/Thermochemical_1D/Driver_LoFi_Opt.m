clear;
close all;
clc;
run('../../src/Set_Paths');

con_hifi = load('HiFi_Opt_Results.mat', 'con').con;
reg_coeff = load('HiFi_Opt_Results.mat', 'obj').obj.reg_coeff;
n_y = con_hifi.fe.m;
n_t = con_hifi.n_t;
n_z = con_hifi.n_z;

con = Thermochemical_LoFi_Constraint_AD(con_hifi);
con.AD_Initialization('Lofi_AD_Files');
obj = Thermochemical_Dynamic_Objective(con, reg_coeff);
opt = Reduced_Space_Optimization(obj, con);

z0 = load('HiFi_Opt_Results.mat', 'z').z;
[u, z] = opt.Optimize(z0);

f = con.con_hifi.Map_Controller_to_Mesh(z);
u_rs = reshape(u, 4 * n_y, n_t);
T = con.I_T * u_rs;
u1 = con.I_u1 * u_rs;
v1 = con.I_v1 * u_rs;
v2 = con.I_v2 * u_rs;

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
surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', v1);
xlabel('Space');
ylabel('Time');
title('v1');
colorbar();
set(gca, 'fontsize', 18);

figure(4);
surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', v2);
xlabel('Space');
ylabel('Time');
title('v2');
colorbar();
set(gca, 'fontsize', 18);

figure(5);
surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', f);
xlabel('Space');
ylabel('Time');
title('f');
colorbar();
set(gca, 'fontsize', 18);

con.Clear_AD();
save('LoFi_Opt_Results.mat', 'obj', 'con', 'u', 'z');
