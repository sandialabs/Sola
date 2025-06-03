clear;
close all;
clc;
addpath(genpath('../../src'));

m = 200;
diff_coeff = 1;
vel_coeff = 1 / 2;
robin_coeff = 2;
reg_coeff = 10;
xi = 1;
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff(m, diff_coeff, vel_coeff, xi);
con_lofi = Diff(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_hifi.x;
z0 = rand(m, 1);

mms_check = true;
grid_refinement_check = true;
finite_diff_check = true;

if mms_check
    y = cos(2 * pi * x);
    z = (10^-2) * (4*pi^2) * cos(2 * pi * x);
    u = con_lofi.State_Solve(z);

    figure;
    hold on;
    plot(x, u, 'LineWidth', 3);
    plot(x, y, '--', 'LineWidth', 3);

    z = (10^-2) * ( (4*pi^2) * cos(2 * pi * x) - (2*pi) * sin(2*pi*x) * vel_coeff );
    u = con_hifi.State_Solve(z);

    figure;
    hold on;
    plot(x, u, 'LineWidth', 3);
    plot(x, y, '--', 'LineWidth', 3);
end

if grid_refinement_check
    m_mesh = 2.^(4:10);
    N = length(m_mesh);
    error = zeros(N, 1);
    for k = 1:N
        con_k = Adv_Diff(m_mesh(k), diff_coeff, vel_coeff, xi);
        x = con_k.x;
        y = cos(2 * pi * x);
        z = (10^-2) * ( (4*pi^2) * cos(2 * pi * x) - (2*pi) * sin(2*pi*x) * vel_coeff );
        u = con_k.State_Solve(z);
        error(k) = sqrt((y - u)' * con_k.M * (y - u));
    end
    figure;
    loglog(1 ./ m_mesh, error);
end

if finite_diff_check
    diffs = opt_lofi.Finite_Difference_Gradient_Check(z0);
    diffs = opt_lofi.Finite_Difference_Hessian_Check(z0);
    diffs = opt_hifi.Finite_Difference_Gradient_Check(z0);
    diffs = opt_hifi.Finite_Difference_Hessian_Check(z0);
end
