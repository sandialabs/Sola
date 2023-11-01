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
[V, D] = eig(linsolve(gevp.A, eye(m)), gevp.M, 'vector');
V = V(:, 1:num_evals);
D = D(1:num_evals);
for k = 1:num_evals
    V(:, k) = sign(evecs(1, k)) * sign(V(1, k)) * V(:, k) / sqrt(V(:, k)' * gevp.M * V(:, k));
end

if ~suppress_figures
    x = linspace(0, 1, m)';
    k = 10;
    figure;
    plot(x, V(:, k), x, evecs(:, k));
end

evecs_ref = load('reference_solution.mat', 'evecs').evecs;
evals_ref = load('reference_solution.mat', 'evals').evals;
ref_diff = max(norm(evals_ref - evals), norm(evecs_ref - evecs));
