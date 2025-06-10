clear;
close all;
clc;
addpath(genpath('../../src'));

diff_coeff = 1;
adv_coeff = 5;
Hmax = .05;

pde_meshing = PDE_Meshing(Hmax);
adv_diff = Adv_Diff(pde_meshing, diff_coeff, adv_coeff);

save('Assembled_Operators.mat');
