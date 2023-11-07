clear;
close all;
clc;
run('../../src/Set_Paths');

m = 50;
n = m;
T = 1;
N = 100;

con = Transient_Inv_Prob_Constraint_AD(m, n, T, N);
forcing = @(x, t) exp(-100 * (x - 0.5).^2);
z_true = (1.e-2) * (2 + cos(3 * pi * con.x));

Generate_Obs_Data(con, z_true, forcing);
