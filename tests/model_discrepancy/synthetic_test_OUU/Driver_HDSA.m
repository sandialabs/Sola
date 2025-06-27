%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

suppress_figures = true;

m = 51;

u_prior_interface = MD_u_Prior_Interface_synthetic_test_OUU(m);
z_prior_interface = MD_z_Prior_Interface_synthetic_test_OUU(m);
