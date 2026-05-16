clear;
close all;

rng(305);

budget = 15;
d = 100;

for i = 1:10
    costs = randi([1, 5], 100, 1);
    vals = costs + randi([0, 5], d, 1);
    f = @(S) sum(vals(S));
    
    eps = 0.000001 ;
    num_trials = 5;
    random_vals = [];
    [~, greedy_val] = Lazy_Greedy_Solve_Knapsack_Cons(f, costs, d, budget);
    
    for j = 1:num_trials-1
        [~, val] = Lazy_Stochastic_CB_Greedy_Solve(f, costs, d, budget, eps);
        random_vals(j) = val;
    end
    
    % run once we know it will break
    [~, val] = Lazy_Stochastic_CB_Greedy_Solve(f, costs, d, budget, eps);
    random_vals(num_trials) = val;
    
    mean_val = mean(random_vals);
    fprintf('Greedy Value=%.4f; Mean Stochastic Value=%.4f\n', greedy_val, mean_val);
end