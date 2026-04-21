%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));

diff_coeff = 2;
adv_coeff = 3;

pde_meshing = PDE_Meshing(.1);
adv_diff = Adv_Diff(pde_meshing, diff_coeff, adv_coeff);
nonlinear_adv_diff = Nonlinear_Adv_Diff(adv_diff);

x = pde_meshing.x;
y = pde_meshing.y;
u_true = (-3 - 2 * x + x.^2) .* (-3 - 2 * y + y.^2);
z = diff_coeff * (-2 * (-3 - 2 * x + x.^2) - 2 * (-3 - 2 * y + y.^2)) + adv_coeff * ((-2 + 2 * x) .* (-3 - 2 * y + y.^2) + (-2 + 2 * y) .* (-3 - 2 * x + x.^2));
u = adv_diff.State_Solve(z);

e = u - u_true;
error = sqrt(e' * pde_meshing.M * e) / sqrt(u_true' * pde_meshing.M * u_true);
disp(['Relative error in state solution = ', num2str(error)]);

plot_fields = true;
if plot_fields
    pde_meshing.Plot_Field(u_true, 'True Solution');
    pde_meshing.Plot_Field(u, 'Low-fidelity Solution');
end

z = diff_coeff * (-2 * (-3 - 2 * x + x.^2) - 2 * (-3 - 2 * y + y.^2)) + adv_coeff * u_true .* ((-2 + 2 * x) .* (-3 - 2 * y + y.^2) + (-2 + 2 * y) .* (-3 - 2 * x + x.^2));

u_hifi = nonlinear_adv_diff.State_Solve(z);
e = u_hifi - u_true;
error = sqrt(e' * pde_meshing.M * e) / sqrt(u_true' * pde_meshing.M * u_true);
disp(['Relative error in state solution = ', num2str(error)]);
if plot_fields
    pde_meshing.Plot_Field(u_hifi, 'High-fidelity Solution');
end
