%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));

n_y = 50;
n_z = n_y;
T = 1;
n_t = 100;

con = Thermal_Constraint(n_y, n_z, T, n_t);
x = con.x;

jacobian_check = true;
mms_check = true;
grid_refinement_check = true;
finite_diff_check = true;

if jacobian_check
    y = randn(n_y, 1);
    z = randn(n_y, 1);
    t = rand;

    con.f_Jacobian_Check(y, z, t);
    con.f_Hessian_Check(y, z, t);
end

if mms_check
    z = 1 + x.^2;
    con.forcing = @(x, t) (10^-3) * (cos(2 * pi * x) + (2 * pi * t) * (2 * x .* sin(2 * pi * x) + 2 * pi * (1 + x.^2) .* cos(2 * pi * x)));
    y = 1 + cos(2 * pi * con.x) * con.t_mesh';
    u = con.State_Solve(z);
    u_reshape = reshape(u, n_y, n_t);
    error = sqrt((y - u_reshape)' * con.M * (y - u_reshape));
    w = ones(n_t, 1);
    w(2:end - 1) = 2;
    this.w = T * w / sum(w);
    error = diag(error)' * w;
    figure(1);
    for k = 1:n_t
        plot(x, u_reshape(:, k), x, y(:, k), '--', 'LineWidth', 3);
        pause(.05);
    end
    disp(['MMS error = ', num2str(error)]);
end

if grid_refinement_check
    m_mesh = 2.^(4:10);
    n_t = length(m_mesh);
    error = zeros(n_t, 1);
    for k = 1:n_t
        con_k = Thermal_Constraint(m_mesh(k), n_z, T, n_t);
        x = con_k.x;
        z = 1 + x.^2;
        con_k.forcing = @(x, t) (10^-3) * (cos(2 * pi * x) + (2 * pi * t) * (2 * x .* sin(2 * pi * x) + 2 * pi * (1 + x.^2) .* cos(2 * pi * x)));
        y = 1 + cos(2 * pi * con_k.x) * con_k.t_mesh';
        u = con_k.State_Solve(z);
        u_reshape = reshape(u, m_mesh(k), n_t);
        error_tmp = sqrt((y - u_reshape)' * con_k.M * (y - u_reshape));
        error(k) = diag(error_tmp)' * con_k.w;
    end
    figure;
    loglog(1 ./ m_mesh, error);
end

if finite_diff_check
    likelihood = Thermal_Likelihood_Model(n_y, n_t);
    prior = Thermal_Prior_Model(con);
    bayes_inv = Bayesian_Inversion(likelihood, prior, con);
    z0 = 10 * rand(n_y, 1);
    diffs = bayes_inv.opt.Finite_Difference_Gradient_Check(z0);
    diffs = bayes_inv.opt.Finite_Difference_Hessian_Check(z0);
end
