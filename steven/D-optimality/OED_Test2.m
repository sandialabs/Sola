%%
%  OED Test
%% Problem Setup
clear;
close all;
N = 100;
x = linspace(0, 1, N)';
m = 2 * ones(N, 1);

% set up forward
cons = Poisson_Constraint(N, 1);
M = cons.M;
S = cons.S;

% set up prior
prior = Poisson_Prior_Model(cons, 5 / 6, 1 / 30);

% set up observation operator
noise_std = 1e-2;
obs_vec = (1:9:100)';
likelihood = Poisson_Likelihood_Model(1 / noise_std, obs_vec, N);
num_sensors = 5;

%%% OED
num_obs = numel(obs_vec);
%% Compute greedy with solution with fast marginal gains approach
linear_oed = Linear_OED_D_Opt(likelihood, prior, cons, num_sensors);
[sensors, val] = linear_oed.Optimize_Design();

fprintf('EIG=%.4f\n', val);

%% Compute greedy solution using log det in data space
B = zeros(num_obs, num_obs);
for i = 1:num_obs
    temp = zeros(num_obs, 1);
    temp(i) = 1.0;
    temp = likelihood.Observation_Operator_Transpose_Apply(temp);
    temp = cons.c_u_Transpose_Inverse_Apply(temp);
    temp = prior.Prior_Covariance_Apply(temp);
    temp = cons.State_Solve(temp);
    B(i, :) = likelihood.Observation_Operator_Apply(temp);
end

G_noise = noise_std^2 * eye(num_obs);

EIG_func(sensors, B, G_noise);

% TODO: Test that the EIG values are the same

greedy_sensors = Lazy_Greedy_Solve_Cardinality_Cons(@(S) EIG_func(S, B, G_noise), ...
                                                    num_obs, ...
                                                    num_sensors);
EIG_func(greedy_sensors, B, G_noise);

% TODO: Compare the solutions found by the different approaches. Should be
%   the same here.

function eig = EIG_func(S, B, G_noise)
    temp = det((B(S, S)) + G_noise(S, S)) / det(G_noise(S, S));
    eig = 0.5 * log(temp);
end
