%% Sample from the prior
%
clear;
close all;
N = 100;
x = linspace(0, 1, N)';
m = 2 * ones(N, 1);

cons = Poisson_Constraint(N, 1);
M = cons.M;
S = cons.S;

prior = Poisson_Prior_Model(cons, 5/6, 1/30);

figure;
hold

num_samples = 100;
samples = zeros(num_samples, N);
for i = 1:num_samples
    z = randn(N, 1);
    samples(i, :) = prior.Square_Root_Apply(z);
    plot(x, samples(i, :));
end
