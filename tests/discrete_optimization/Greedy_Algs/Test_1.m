clear;
close all;

rng(200);

budget = 5;
num_facilities = 50;
num_customers = 20;
M = rand(num_customers, num_facilities);
costs = ones(num_facilities, 1);

f = @(S) sum(max(M(:, S)));

% Brute force find the optimal solution
combinations = nchoosek(1:num_facilities, budget);
f_vals = zeros(size(combinations, 1), 1);

for i = 1:size(combinations, 1)
    f_vals(i) = f(combinations(i, :));
end
[best_f, best_idx] = max(f_vals);

%% Test Deterministic Greedy
[greedy_sol, greedy_val] = Lazy_Greedy_Solve_Knapsack_Cons(f, costs, num_facilities, budget);

passed = (greedy_val >= (1 - exp(-1)) * best_f);
%% Test Stochastic Greedy
num_trials = 100;
random_vals = zeros(num_trials, 1);
eps = 0.1;

for i = 1:num_trials
    [~, val] = Lazy_Stochastic_CB_Greedy_Solve(f, costs, num_facilities, budget, eps);
    random_vals(i) = val;
end

passed = passed & (mean(random_vals) >= (1 - exp(-1) - eps) * best_f);

if passed
    fprintf(1, '\ndiscrete_optimization/Greedy_Algs passed.\n');
else
    fprintf(2, '\ndiscrete_optimization/Greedy_Algs failed.\n');
end
