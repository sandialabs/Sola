clear
close all
clc
addpath(genpath('../../src'))

m = 50;
con = Darcy_Constraint(m);
forcing = @(x) 1 + 0*x;

Generate_Obs_Data(con,forcing)
