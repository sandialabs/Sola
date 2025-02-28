clear;
close all;
clc;
addpath(genpath('../../src'));

diff_coeff = 2;
adv_coeff = 3;

pde_meshing = PDE_Meshing(.05);
diff = Diff_Constraint(pde_meshing, diff_coeff);
adv_diff = Adv_Diff(pde_meshing, diff_coeff, adv_coeff);

x = pde_meshing.x;
y = pde_meshing.y;
u_true = cos(2 * pi * x) .* (- 2 * y + y.^2);
z = diff_coeff * (4 * pi^2 * cos(2 * pi * x) .* (-2 * y + y.^2) - 2 * cos(2 * pi * x));
z = z + adv_coeff * (-2 * pi * sin(2 * pi * x) .* (-2 * y + y.^2) + (-2 + 2 * y) .* cos(2 * pi * x));
u = adv_diff.State_Solve(z);

e = u - u_true;
error = sqrt(e' * pde_meshing.M * e) / sqrt(u_true' * pde_meshing.M * u_true);
disp(['Relative error in state solution = ', num2str(error)]);

plot_fields = true;
if plot_fields
    pde_meshing.Plot_Field(u_true, 'True Solution');
    pde_meshing.Plot_Field(u, 'Adv-Diff Model Solution');
end

z = diff_coeff * (4 * pi^2 * cos(2 * pi * x) .* (-2 * y + y.^2) - 2 * cos(2 * pi * x));

u_diff = diff.State_Solve(z);
e = u_diff - u_true;
error = sqrt(e' * pde_meshing.M * e) / sqrt(u_true' * pde_meshing.M * u_true);
disp(['Relative error in state solution = ', num2str(error)]);
if plot_fields
    pde_meshing.Plot_Field(u_diff, 'Diff Model Solution');
end
