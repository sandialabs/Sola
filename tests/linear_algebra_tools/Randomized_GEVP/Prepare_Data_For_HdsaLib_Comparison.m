%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
rng(1431242);

random_numbers = randn(10^5, 1);
writematrix(random_numbers, 'random_numbers.txt');

rng(1431242);

suppress_figures = true;

m = 50;
gevp = Randomized_GEVP_Test(m);

num_evals = 20;
oversampling = 20;
[evecs, evals] = gevp.Compute_GEVP(num_evals, oversampling);

save('Sabl_Output.mat', 'evecs', 'evals');
