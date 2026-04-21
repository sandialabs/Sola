%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));

m = 200;
diff_coeff = 1;
react_coeff = -10;
reg_coeff = 1.e-5;
obj = Diff_React_Objective(m, reg_coeff);
con = Diff_React_Constraint(m, diff_coeff, react_coeff);

opt = Reduced_Space_Optimization(obj, con);
x = con.x;
z0 = rand(m, 1);

% We have to be careful with the tests. The PDE admits multiple solutions
% for certain reaction functions. Skipping to a different solution messes
% up the tests

mms_check = true;
grid_refinement_check = true;
finite_diff_check = true;
reaction_function_check = true;

if mms_check
    y = -2 * x.^3 + 3 * x.^2;
    z = -diff_coeff * (-12 * x + 6) - react_coeff * con.Reaction_Function(y);
    u = con.State_Solve(z);

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
        con_k = Diff_React_Constraint(m_mesh(k), diff_coeff, react_coeff);
        x = con_k.x;
        y = -2 * x.^3 + 3 * x.^2;
        z = -diff_coeff * (-12 * x + 6) - react_coeff * con.Reaction_Function(y);
        u = con_k.State_Solve(z);
        error(k) = sqrt((y - u)' * con_k.M * (y - u));
    end
    figure;
    loglog(1 ./ m_mesh, error);
end

if finite_diff_check
    diffs = opt.Finite_Difference_Gradient_Check(z0);
    diffs = opt.Finite_Difference_Hessian_Check(z0);
end

if reaction_function_check
    u = randn(m, 1);
    lambda = randn(m, 1);
    z = randn(m, 1);
    con.Finite_Difference_Reaction_Function_Jacobian(u);
    con.Finite_Difference_Reaction_Function_Hessian(u, lambda);
    con.Finite_Difference_Constraint_Hessian(u, z, lambda);
end
