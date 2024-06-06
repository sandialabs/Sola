clear;
close all;
addpath(genpath('../../../src'));
rng(1431242);

suppress_figures = true;

m = 50;
gevp = Randomized_GEVP_Test(m);

num_evals = 20;
oversampling = 20;
[evecs, evals] = gevp.Compute_GEVP(num_evals, oversampling);

save('Sabl_Output.mat', 'evecs','evals');
