%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

m = 50;
diff_coeff = 1;
con = Poisson_Constraint(m, diff_coeff);
x = con.x;

likelihood = Poisson_Likelihood_Model(con);
prior = Poisson_Prior_Model(con);
bayes_inv = Bayesian_Inversion(likelihood, prior, con);

mms_check = true;
grid_refinement_check = true;
finite_diff_check = true;

if mms_check
    y = x .* (1 - x);
    z = 2 * diff_coeff * ones(m, 1);
    u = con.State_Solve(z);

    figure;
    hold on;
    plot(x, u, 'LineWidth', 3);
    plot(x, y, '--', 'LineWidth', 3);
end

if grid_refinement_check
    m_mesh = 2.^(2:10);
    N = length(m_mesh);
    error = zeros(N, 1);
    for k = 1:N
        con_k = Poisson_Constraint(m_mesh(k), diff_coeff);
        x = con_k.x;
        y = x .* (1 - x);
        z = 2 * diff_coeff * ones(m_mesh(k), 1);
        u = con_k.State_Solve(z);
        error(k) = sqrt((y - u)' * con_k.M * (y - u));
    end
    figure;
    loglog(1 ./ m_mesh, error);
end

if finite_diff_check
    z0 = 100 * randn(m, 1);
    diffs = bayes_inv.opt.Finite_Difference_Gradient_Check(z0);
    diffs = bayes_inv.opt.Finite_Difference_Hessian_Check(z0);
end
