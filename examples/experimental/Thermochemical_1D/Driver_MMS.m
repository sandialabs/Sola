clear;
close all;
clc;

n_y = 50;
T = .2;
n_t = 11;
n_z = n_y * n_t * 6;

con = Thermochemical_Constraint_AD_MMS(n_y, n_z, T, n_t);
con.AD_Initialization();

%% MMS test
x = con.fe.x;
t = con.t_mesh;

T = 1 + cos(2 * pi * x) * t';
u1 = cos(4 * pi * x) * t';
u2 = cos(2 * pi * x) * t';
v1 = 20 + cos(2 * pi * x) * t';
v2 = 10 + cos(4 * pi * x) * t';
v3 = 12 + cos(2 * pi * x) * t';

T_dot = cos(2 * pi * x) * ones(1, n_t);
u1_dot = cos(4 * pi * x) * ones(1, n_t);
u2_dot = cos(2 * pi * x) * ones(1, n_t);
v1_dot = cos(2 * pi * x) * ones(1, n_t);
v2_dot = cos(4 * pi * x) * ones(1, n_t);
v3_dot = cos(2 * pi * x) * ones(1, n_t);

T_lap = -con.diff_T * 4 * pi^2 * cos(2 * pi * x) * t';
u1_lap = -con.diff_u1 * 16 * pi^2 * cos(4 * pi * x) * t';
u2_lap = -con.diff_u2 * 4 * pi^2 * cos(2 * pi * x) * t';
v1_lap = -con.diff_v1 * 4 * pi^2 * cos(2 * pi * x) * t';
v2_lap = -con.diff_v2 * 16 * pi^2 * cos(4 * pi * x) * t';
v3_lap = -con.diff_v3 * 4 * pi^2 * cos(2 * pi * x) * t';

R = con.Evaluate_Arrhenius_Law(T);

f_T = T_dot - T_lap + con.cooling_1 * con.react_rate_1 * (R .* v1 .* v2) + con.cooling_2 * con.react_rate_2 * (R .* v1 .* v3);
f_u1 = u1_dot - u1_lap - con.react_rate_1 * (R .* v1 .* v2);
f_u2 = u2_dot - u2_lap - con.react_rate_2 * (R .* v1 .* v3);
f_v1 = v1_dot - v1_lap + con.react_rate_1 * (R .* v1 .* v2) + con.react_rate_2 * (R .* v1 .* v3);
f_v2 = v2_dot - v2_lap + con.react_rate_1 * (R .* v1 .* v2);
f_v3 = v3_dot - v3_lap + con.react_rate_2 * (R .* v1 .* v3);

tmp = [f_T; f_u1; f_u2; f_v1; f_v2; f_v3];
z = tmp(:);

u = con.State_Solve(z);

u_rs = reshape(u, 6 * n_y, n_t);
T_sol = u_rs(1:n_y, :);
u1_sol = u_rs((n_y + 1):(2 * n_y), :);
u2_sol = u_rs((2 * n_y + 1):(3 * n_y), :);
v1_sol = u_rs((3 * n_y + 1):(4 * n_y), :);
v2_sol = u_rs((4 * n_y + 1):(5 * n_y), :);
v3_sol = u_rs((5 * n_y + 1):(6 * n_y), :);

figure(1);
for k = 1:n_t
    plot(x, T(:, k), x, T_sol(:, k), '--', 'LineWidth', 3);
    pause(0.05);
end

figure(2);
for k = 1:n_t
    plot(x, u1(:, k), x, u1_sol(:, k), '--', 'LineWidth', 3);
    pause(0.05);
end

figure(3);
for k = 1:n_t
    plot(x, u2(:, k), x, u2_sol(:, k), '--', 'LineWidth', 3);
    pause(0.05);
end

figure(4);
for k = 1:n_t
    plot(x, v1(:, k), x, v1_sol(:, k), '--', 'LineWidth', 3);
    pause(0.05);
end

figure(5);
for k = 1:n_t
    plot(x, v2(:, k), x, v2_sol(:, k), '--', 'LineWidth', 3);
    pause(0.05);
end

figure(6);
for k = 1:n_t
    plot(x, v3(:, k), x, v3_sol(:, k), '--', 'LineWidth', 3);
    pause(0.05);
end
