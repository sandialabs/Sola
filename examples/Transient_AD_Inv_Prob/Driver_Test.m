clear;
close all;
clc;
run('../../src/Set_Paths');

m = 50;
n = m;
T = 1;
N = 100;

con = Transient_Inv_Prob_Constraint_AD(m, n, T, N);
x = con.x;

jacobian_check = false;
mms_check = false;
grid_refinement_check = false;
finite_diff_check = true;

if jacobian_check
    y = randn(m, 1);
    z = randn(m, 1);
    t = rand;

    con.AD_Initialization();
    diffs = con.Time_Instance_RHS_Jacobian_y_Check(y, z, t);
    diffs = con.Time_Instance_RHS_Jacobian_z_Check(y, z, t);
    diffs = con.Time_Instance_RHS_Hessian_yz_Check(y, z, t);
    diffs = con.Time_Instance_RHS_Hessian_zy_Check(y, z, t);
end

if mms_check
    z = 1 + x.^2;
    con.forcing = @(x, t) (10^-3) * (cos(2 * pi * x) + (2 * pi * t) * (2 * x .* sin(2 * pi * x) + 2 * pi * (1 + x.^2) .* cos(2 * pi * x)));
    con.AD_Initialization();
    y = 1 + cos(2 * pi * con.x) * con.t_mesh';
    u = con.State_Solve(z);
    u_reshape = reshape(u, m, N);
    error = sqrt((y - u_reshape)' * con.M * (y - u_reshape));
    error = diag(error)' * con.w;
    figure(1);
    for k = 1:N
        plot(x, u_reshape(:, k), x, y(:, k), '--', 'LineWidth', 3);
        pause(.05);
    end
    disp(['MMS error = ', num2str(error)]);
end

if grid_refinement_check
    m_mesh = 2.^(4:10);
    N = length(m_mesh);
    error = zeros(N, 1);
    for k = 1:N
        con_k = Transient_Inv_Prob_Constraint_AD(m_mesh(k), m_mesh(k), T, N);
        x = con_k.x;
        z = 1 + x.^2;
        con_k.forcing = @(x, t) (10^-3) * (cos(2 * pi * x) + (2 * pi * t) * (2 * x .* sin(2 * pi * x) + 2 * pi * (1 + x.^2) .* cos(2 * pi * x)));
        con_k.AD_Initialization();
        y = 1 + cos(2 * pi * con_k.x) * con_k.t_mesh';
        u = con_k.State_Solve(z);
        u_reshape = reshape(u, m_mesh(k), N);
        error_tmp = sqrt((y - u_reshape)' * con_k.M * (y - u_reshape));
        error(k) = diag(error_tmp)' * con_k.w;
    end
    figure;
    loglog(1 ./ m_mesh, error);
end

if finite_diff_check
    con.AD_Initialization();
    likelihood = Transient_Inv_Prob_Likelihood_Model(m, N);
    prior = Transient_Inv_Prob_Prior_Model(con);
    bayes_inv = Bayesian_Inversion(likelihood, prior, con);
    z0 = 10 * rand(m, 1);
    diffs = bayes_inv.opt.Finite_Difference_Gradient_Check(z0);
    diffs = bayes_inv.opt.Finite_Difference_Hessian_Check(z0);
end
