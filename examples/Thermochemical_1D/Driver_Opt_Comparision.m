clear;
close all;
clc;

reg_coeff = load('HiFi_Opt_Results.mat', 'obj').obj.reg_coeff;

con_hifi = load('HiFi_Opt_Results.mat', 'con').con;
con_hifi.AD_Initialization('Hifi_AD_Files');
obj_hifi = Thermochemical_Dynamic_Objective(con_hifi, reg_coeff);
opt_hifi = Reduced_Space_Optimization(obj_hifi, con_hifi);
u_hifi = load('HiFi_Opt_Results.mat', 'u').u;
z_hifi = load('HiFi_Opt_Results.mat', 'z').z;

u_lofi = load('LoFi_Opt_Results.mat', 'u').u;
z_lofi = load('LoFi_Opt_Results.mat', 'z').z;

val_hifi = opt_hifi.Jhat(z_hifi);
val_lofi = opt_hifi.Jhat(z_lofi);
state_hifi = u_hifi;
state_lofi = con_hifi.State_Solve(z_lofi);

f_hifi = con_hifi.Map_Controller_to_Mesh(z_hifi);
f_lofi = con_hifi.Map_Controller_to_Mesh(z_lofi);

figure(1);
surf(con_hifi.fe.x * ones(1, con_hifi.n_t), ones(con_hifi.fe.m, 1) * con_hifi.t_mesh', f_hifi);
xlabel('Space');
ylabel('Time');
title('HiFi f');
colorbar();
set(gca, 'fontsize', 18);

figure(2);
surf(con_hifi.fe.x * ones(1, con_hifi.n_t), ones(con_hifi.fe.m, 1) * con_hifi.t_mesh', f_lofi);
xlabel('Space');
ylabel('Time');
title('LoFi f');
colorbar();
set(gca, 'fontsize', 18);

con_hifi.Clear_AD();
