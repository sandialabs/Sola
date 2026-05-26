%%
%  OED Test
%% Problem Setup
clear;
close all;
N = 100;    % parameter dimension
x = linspace(0, 1, N)';
m = 2 * ones(N, 1);

% set up forward
cons = Poisson_Constraint(N, 1);
M = cons.M;
S = cons.S;

% set up prior
prior = Poisson_Prior_Model(cons, 5 / 6, 1 / 30);

% set up observation model
noise_std = 1e-2;
obs_vec = (1:9:100)';
likelihood = Poisson_Likelihood_Model(noise_std, obs_vec, N);
num_sensors = 5;

%% OED
num_obs = numel(obs_vec);
%% Compute greedy with solution with data space log det approach
linear_oed = Linear_OED_D_Opt(likelihood, prior, cons, num_sensors);
greedy_sensors = linear_oed.Optimize_Design();

%% Compute greedy solution using fast marginal gain computations
% construct F
F = zeros(num_obs, N);
for i = 1:num_obs
    x = zeros(num_obs, 1);
    x(i) = 1;
    temp = likelihood.Observation_Operator_Transpose_Apply(x);
    temp = cons.c_u_Transpose_Inverse_Apply(temp) / noise_std;
    temp = -cons.c_z_Transpose_Apply(temp);
    temp = prior.Mass_Matrix_Inverse_Apply(temp);
    f = prior.Prior_Covariance_Factor_Apply(temp);
    F(i, :) = f;
end

sensors = zeros(1, num_sensors);
C = RankOneUpdatesMatrix(M);
val = 0;

% initialize list of marginal gains
marginal_gains = RedBlackTree();
for v = 1:num_obs
    marginal_gains.Insert(inf, ...
                          struct('sensor', v, 'time', 0));
end

for i = 1:num_sensors
    while true
        [gain, data] = marginal_gains.PopMax();
        v = data.sensor;
        t = data.time;
        f = F(v, :)';
        if t == i
            sensors(i) = v;
            C.Add_Update(f);
            val = val + gain;
            break
        else
            g = C.Inverse_Apply(f);
            gain = log(1 + g' * M * f);
            marginal_gains.Insert(gain, ...
                                  struct('sensor', v, 'time', i));
        end
    end
end

fast_mg_eig = 0.5 * val;
data_space_eig = linear_oed.OED_Objective(greedy_sensors);

error = abs(fast_mg_eig - data_space_eig);

if error > 1.e-9
    fprintf(2, '\noptimal_experimental_design/Poisson_D_Opt failed.\n');
else
    fprintf(1, '\noptimal_experimental_design/Poisson_D_Opt passed.\n');
end
