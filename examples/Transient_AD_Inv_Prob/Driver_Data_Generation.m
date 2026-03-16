clear;
close all;
clc;

n_y = 50;
n_z = n_y;
T = 1;
n_t = 100;

con = Transient_Inv_Prob_Constraint_AD(n_y, n_z, T, n_t);
forcing = @(x, t) exp(-100 * (x - 0.5).^2);
z_true = (1.e-2) * (2 + cos(3 * pi * con.x));

Generate_Obs_Data(con, z_true, forcing);
