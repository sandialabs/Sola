%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));

m = 200;
mu = 1;
gamma = 0.5;
reg_coeff = 1.e-4;
obj = Heat_Problem_Objective(m, reg_coeff);
con_lofi = Heat_Problem_LoFi_Constraint(m, mu);
con_hifi = Heat_Problem_HiFi_Constraint(con_lofi, gamma);

opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;
z0 = rand(m, 1);

mms_check = true;
grid_refinement_check = true;
finite_diff_check = true;
reaction_function_check = false;
diffusion_function_check = false;

% Manufactured solution y = sin(pi x), which satisfies y(0) = y(1) = 0.
mms_solution = @(xx) sin(pi * xx);
mms_lofi_source = @(xx) pi^2 * sin(pi * xx) + mu * sin(pi * xx).^3;
mms_hifi_source = @(xx) pi^2 * sin(pi * xx) ...
    - gamma * pi^2 * (2 * sin(pi * xx) .* cos(pi * xx).^2 - sin(pi * xx).^3) ...
    + mu * sin(pi * xx).^3;

if mms_check
    y = mms_solution(x);

    z = mms_lofi_source(x);
    u = con_lofi.State_Solve(z);
    figure;
    hold on;
    plot(x, u, 'LineWidth', 3);
    plot(x, y, '--', 'LineWidth', 3);
    title('Low-fidelity MMS', 'Interpreter', 'latex');
    legend({'$u$', 'exact'}, 'Interpreter', 'latex');

    z = mms_hifi_source(x);
    u = con_hifi.State_Solve(z);
    figure;
    hold on;
    plot(x, u, 'LineWidth', 3);
    plot(x, y, '--', 'LineWidth', 3);
    title('High-fidelity MMS', 'Interpreter', 'latex');
    legend({'$u$', 'exact'}, 'Interpreter', 'latex');
end

if grid_refinement_check
    m_mesh = 2.^(4:10);
    N = length(m_mesh);
    error_lofi = zeros(N, 1);
    error_hifi = zeros(N, 1);
    for k = 1:N
        con_k_lofi = Heat_Problem_LoFi_Constraint(m_mesh(k), mu);
        con_k_hifi = Heat_Problem_HiFi_Constraint(con_k_lofi, gamma);
        xk = con_k_lofi.x;
        y = mms_solution(xk);

        u = con_k_lofi.State_Solve(mms_lofi_source(xk));
        error_lofi(k) = sqrt((y - u)' * con_k_lofi.M * (y - u));

        u = con_k_hifi.State_Solve(mms_hifi_source(xk));
        error_hifi(k) = sqrt((y - u)' * con_k_hifi.M * (y - u));
    end
    figure;
    loglog(1 ./ m_mesh, error_lofi, '-o', 1 ./ m_mesh, error_hifi, '-s');
    xlabel('$h$', 'Interpreter', 'latex');
    ylabel('$L^2$ error', 'Interpreter', 'latex');
    legend({'low-fidelity', 'high-fidelity'}, 'Interpreter', 'latex');
end

if finite_diff_check
    diffs = opt_lofi.Finite_Difference_Gradient_Check(z0);
    diffs = opt_lofi.Finite_Difference_Hessian_Check(z0);
    diffs = opt_hifi.Finite_Difference_Gradient_Check(z0);
    diffs = opt_hifi.Finite_Difference_Hessian_Check(z0);
end

if reaction_function_check
    u = randn(m, 1);
    lambda = randn(m, 1);
    con_lofi.Finite_Difference_Reaction_Function_Jacobian(u);
    con_lofi.Finite_Difference_Reaction_Function_Hessian(u, lambda);
    con_hifi.Finite_Difference_Reaction_Function_Jacobian(u);
    con_hifi.Finite_Difference_Reaction_Function_Hessian(u, lambda);
end

if diffusion_function_check
    u = randn(m, 1);
    lambda = randn(m, 1);
    z = randn(m, 1);
    con_hifi.Finite_Difference_Diffusion_Function_Jacobian(u);
    con_hifi.Finite_Difference_Diffusion_Function_Hessian(u, lambda);
    con_hifi.Finite_Difference_Constraint_Hessian(u, z, lambda);
end
