clear;
close all;

rng(200);
vals = randi([1, 10], 100, 1);
f = @(S) sum(vals(S));

budget = 10;
d = 100;
costs = randi([1, 10], d, 1);
num_trials = 100;
random_vals = [];
[greedy_sol, greedy_val] = Lazy_Greedy_Solve_Knapsack_Cons(f, costs, d, budget);

fprintf('Greedy Value=%.4f;\n', greedy_val);