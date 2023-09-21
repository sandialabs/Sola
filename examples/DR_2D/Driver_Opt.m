clear
close all
clc
addpath(genpath('../../src'))

h = .1;
mesh = PDE_Meshing(h);
x = mesh.x;
y = mesh.y;
M = mesh.M;
diff_react_lofi = Diff_React_Lofi(mesh);
diff_react_hifi = Diff_React_Hifi(mesh);
m = length(x);

reg_coeff = 1.e-4;
diff_react_opt = Diff_React_Opt(diff_react_lofi,reg_coeff);

z0 = rand(m,1);
[u_lofi,z_lofi] = diff_react_opt.Optimize(z0);

T = diff_react_opt.T;
u = diff_react_hifi.State_Solve(diff_react_hifi.Map_z_to_Control_Fun(z_lofi));

name = 'High-fidelity state';
mesh.Plot_Field(u,name)

name = 'Low-fidelity state';
mesh.Plot_Field(u_lofi,name)

name = 'Target state';
mesh.Plot_Field(T,name)

name = 'Low-fidelity control';
mesh.Plot_Field(z_lofi,name)

Z = zeros(length(z_lofi),2);
Y = zeros(length(u),2);

Z(:,1) = z_lofi;
Y(:,1) = u-u_lofi;

Z(:,2) = 0.5*z_lofi + 0.5*250*exp(-10*(x.^2 + y.^2));
Y(:,2) = diff_react_hifi.State_Solve(diff_react_hifi.Map_z_to_Control_Fun(Z(:,2))) - diff_react_lofi.State_Solve(Z(:,2));

save('Optimization_Results.mat','h','reg_coeff','z_lofi','u_lofi','Z','Y')