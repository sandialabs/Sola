clear
close all
clc
addpath(genpath('../../src'))

m = 50;
con = Thermal_Constraint(m);
x = con.x;

jacobian_check = true;
mms_check = true;
grid_refinement_check = true;
finite_diff_check = true;

if jacobian_check
   u = randn(m,1);
   z = randn(m,1);
   [diffs_z,jacobian_z_transpose_check,diffs_u,jacobian_u_transpose_check,solve_res] = con.Finite_Difference_Constraint_Check(u,z);
end

if mms_check
    z = x.^2;
    con.forcing = (10^-3)*(6*x.^2 - 2*x);
    y = x-x.^2;
    u = con.State_Solve(z);
    error = sqrt((y-u)'*con.M*(y-u));
    figure,
    hold on
    plot(x,u,'LineWidth',3)
    plot(x,y,'--','LineWidth',3)
    disp(['MMS error = ',num2str(error)])
end

if grid_refinement_check
    m_mesh = 2.^(4:10);
    N = length(m_mesh);
    error = zeros(N,1);
    for k = 1:N
        con_k = Thermal_Constraint(m_mesh(k));
        x = con_k.x;
        z = x.^2;
        con_k.forcing = (10^-3)*(6*x.^2 - 2*x);
        y = x-x.^2;
        u = con_k.State_Solve(z);
        error(k) = sqrt((y-u)'*con_k.M*(y-u));
    end
    figure,
    loglog(1./m_mesh,error)
end

if finite_diff_check
    likelihood = Thermal_Likelihood_Model(m);
    prior = Thermal_Prior_Model(con);
    bayes_inv = Bayesian_Inversion(likelihood,prior,con);
    z0 = randn(m,1);
    diffs = bayes_inv.opt.Finite_Difference_Gradient_Check(z0);
    diffs = bayes_inv.opt.Finite_Difference_Hessian_Check(z0);
end

