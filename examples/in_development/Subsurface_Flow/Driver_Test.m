%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));

m = 200;
k0 = 1;
alpha = 0.5;
p0 = 0;
viscosity = 1;
reg_coeff = 1.e-6;
hifi_leakoff_coeff = 30;
hifi_leakoff_center = 0.55;
hifi_leakoff_width = 0.25;

obj = Subsurface_Objective(m, reg_coeff, p0);
con_lofi = Subsurface_LoFi_Constraint(m, k0, alpha, p0, viscosity);
con_hifi = Subsurface_HiFi_Constraint(con_lofi, hifi_leakoff_coeff, hifi_leakoff_center, hifi_leakoff_width);

opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;
z0 = 0.1 * randn(m, 1);

mms_check = true;
grid_refinement_check = true;
finite_diff_check = true;
diffusion_function_check = true;

% Manufactured pressure p = sin(pi x), satisfying p(0)=p(1)=p0 when p0=0.
mms_solution = @(xx) p0 + sin(pi * xx);
mms_lofi_source = @(xx) Subsurface_MMS_Source(xx, k0, alpha, p0, viscosity, 'lofi', hifi_leakoff_coeff, hifi_leakoff_center, hifi_leakoff_width);
mms_hifi_source = @(xx) Subsurface_MMS_Source(xx, k0, alpha, p0, viscosity, 'hifi', hifi_leakoff_coeff, hifi_leakoff_center, hifi_leakoff_width);

if mms_check
    p_exact = mms_solution(x);

    z = mms_lofi_source(x);
    u = con_lofi.State_Solve(z);
    figure;
    hold on;
    plot(x, u, 'LineWidth', 3);
    plot(x, p_exact, '--', 'LineWidth', 3);
    title('Low-fidelity MMS', 'Interpreter', 'latex');
    legend({'$p$', 'exact'}, 'Interpreter', 'latex');

    z = mms_hifi_source(x);
    u = con_hifi.State_Solve(z);
    figure;
    hold on;
    plot(x, u, 'LineWidth', 3);
    plot(x, p_exact, '--', 'LineWidth', 3);
    title('High-fidelity MMS', 'Interpreter', 'latex');
    legend({'$p$', 'exact'}, 'Interpreter', 'latex');
end

if grid_refinement_check
    m_mesh = 2.^(4:10);
    N = length(m_mesh);
    error_lofi = zeros(N, 1);
    error_hifi = zeros(N, 1);
    for k = 1:N
        con_k_lofi = Subsurface_LoFi_Constraint(m_mesh(k), k0, alpha, p0, viscosity);
        con_k_hifi = Subsurface_HiFi_Constraint(con_k_lofi, hifi_leakoff_coeff, hifi_leakoff_center, hifi_leakoff_width);
        xk = con_k_lofi.x;
        p_exact = mms_solution(xk);

        u = con_k_lofi.State_Solve(mms_lofi_source(xk));
        error_lofi(k) = sqrt((p_exact - u)' * con_k_lofi.M * (p_exact - u));

        u = con_k_hifi.State_Solve(mms_hifi_source(xk));
        error_hifi(k) = sqrt((p_exact - u)' * con_k_hifi.M * (p_exact - u));
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

if diffusion_function_check
    u = p0 + 0.25 * randn(m, 1);
    u([1, end]) = p0;
    lambda = randn(m, 1);
    z = randn(m, 1);
    con_lofi.Finite_Difference_Diffusion_Function_Jacobian(u);
    con_lofi.Finite_Difference_Diffusion_Function_Hessian(u, lambda);
    con_lofi.Finite_Difference_Constraint_Hessian(u, z, lambda);
    con_hifi.Finite_Difference_Diffusion_Function_Jacobian(u);
    con_hifi.Finite_Difference_Diffusion_Function_Hessian(u, lambda);
    con_hifi.Finite_Difference_Leakoff_Function_Jacobian(u);
    con_hifi.Finite_Difference_Leakoff_Function_Hessian(u, lambda);
    con_hifi.Finite_Difference_Constraint_Hessian(u, z, lambda);
end

function [q] = Subsurface_MMS_Source(x, k0, alpha, p0, viscosity, model, hifi_leakoff_coeff, hifi_leakoff_center, hifi_leakoff_width)
    p = p0 + sin(pi * x);
    px = pi * cos(pi * x);
    pxx = -pi^2 * sin(pi * x);
    switch model
        case 'hifi'
            r = p - p0;
            s = exp(-((x - hifi_leakoff_center) / hifi_leakoff_width).^2);
            a = (k0 / viscosity) * exp(alpha * r);
            a_prime = alpha * a;
            a_x = zeros(size(x));
            leakoff = hifi_leakoff_coeff * s .* (r.^3);
        case 'lofi'
            r = p - p0;
            a = (k0 / viscosity) * exp(alpha * r);
            a_prime = alpha * a;
            a_x = zeros(size(x));
            leakoff = zeros(size(x));
        otherwise
            error('Unknown model.');
    end
    q = -(a_prime .* (px.^2) + a_x .* px + a .* pxx) + leakoff;
end
