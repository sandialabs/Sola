clear;
close all;
clc;
addpath(genpath('../../src'));

m = 50;
n = m;
T = 10;
N = 100;

con = Thermal_Constraint(m, n, T, N);
forcing = @(x, t) exp(-100 * (x - 0.5).^2);
z_true = (1.e-2) * (2 + cos(3 * pi * con.x));

Generate_Obs_Data(con, z_true, forcing);
