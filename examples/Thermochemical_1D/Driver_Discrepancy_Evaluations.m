clear;
close all;
clc;
run('../../src/Set_Paths');

z_lofi = load('LoFi_Opt_Results.mat', 'z').z;

Z = zeros(length(z_lofi), 1);
Z(:, 1) = z_lofi;

con_hifi = load('HiFi_Opt_Results.mat', 'con').con;
con_hifi.AD_Initialization('Hifi_AD_Files');

u_hifi = zeros(con_hifi.n_y * con_hifi.n_t, size(Z, 2));
for k = 1:size(Z, 2)
    u_hifi(:, k) = con_hifi.State_Solve(Z(:, k));
end

con_hifi.Clear_AD();

con = Thermochemical_LoFi_Constraint_AD(con_hifi);
con.AD_Initialization('Lofi_AD_Files');

u_lofi = zeros(con.n_y * con.n_t, size(Z, 2));
for k = 1:size(Z, 2)
    u_lofi(:, k) = con.State_Solve(Z(:, k));
end

D = 0 * u_lofi;
for k = 1:size(D, 2)
    u_hifi_rs = reshape(u_hifi(:, k), 6 * con.fe.m, con.n_t);
    u_hifi_rs = u_hifi_rs(kron([1; 2; 4; 5], (1:con.fe.m)'), :);
    D(:, k) = u_hifi_rs(:) - u_lofi(:, k);
end

save('Discrepancy_Evaluations.mat', 'Z', 'D');

con.Clear_AD();
