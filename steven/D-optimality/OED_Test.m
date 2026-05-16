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
prior = Poisson_Prior_Model(cons, 5/6, 1/30);

% set up observation operator
noise_std = 1e-2;
obs_vec = (1:9:100)';
likelihood = Poisson_Likelihood_Model(noise_std, obs_vec, N);

%% OED
% First: construct the rows of the matrix -- this isn't necessarily what I
% want do. But this is just an initial step
num_obs = numel(obs_vec);
rows = zeros(num_obs, N);

for i = 1:num_obs
    x = zeros(num_obs, 1);
    x(i) = 1;
    temp = likelihood.Observation_Operator_Transpose_Apply(x);
    temp = cons.c_u_Transpose_Inverse_Apply(temp) / noise_std;
    f = prior.Prior_Covariance_Factor_Apply(temp);
    rows(i, :) = f;
end

% Second: what do we call this matrix?
sensors = [];
k = 5;
C = RankOneUpdatesMatrix(M);
val = 0;

% initialize list of marginal gains
marginal_gains = RedBlackTree();
for v = 1:num_obs
    marginal_gains.Insert(inf, ...
        struct('sensor', v, 'time', 0));
end

for i = 1:k
    while true
        [gain, data] = marginal_gains.PopMax();
        v = data.sensor;
        t = data.time;
        f = rows(v, :)';
        if t == i
            sensors(i) = v;
            C.Add_Update(f);
            val = val + gain;
            break;
        else
            g = C.Inverse_Apply(f);
            gain = log(1 + g' * M * f);
            marginal_gains.Insert(gain, ...
                struct('sensor', v, 'time', i)); 
        end
    end
end

fprintf('EIG=%.4f\n', 0.5 * val);
%% plot the sensors
% figure;
% obs_x = obs_vec / 100;
% xline(obs_x(sensors))

%% Construct matrix for log det in data space
B = zeros(num_obs, num_obs);
for i = 1:num_obs
    f = noise_std * rows(i, :)';
    f = prior.Prior_Covariance_Factor_Apply(f);
    f = cons.State_Solve(f);
    B(i, :) = likelihood.Observation_Operator_Apply(f);
end

G_noise = noise_std^2 * eye(num_obs);

function eig = EIG_func(S, B, G_noise)
    temp = det((B(S, S)) + G_noise(S, S)) / det(G_noise(S, S));
    eig = 0.5 * log(temp);
end

EIG_func(sensors, B, G_noise)

greedy_sensors = Lazy_Greedy_Solve_Cardinality_Cons(@(S) EIG_func(S, B, G_noise), ...
    num_obs, k);
EIG_func(greedy_sensors, B, G_noise)