clear;
close all;
rng(142);

print_output = false;

m = 50;
n = m;
T = 1;
N = 10;

obj = Thermal_Objective(m, n, T, N);
obj_AD = Thermal_Objective_AD(m, n, T, N);
obj_AD.verbose = print_output;
evalc('obj_AD.AD_Initialization()');

con = Thermal_Constraint(m, n, T, N);
con_AD = Thermal_Constraint_AD(m, n, T, N);
con_AD.verbose = print_output;
evalc('con_AD.AD_Initialization()');

y = randn(m, 1);
z = randn(n, 1);
t = rand;
lambda = randn(m, 1);
v = randn(m, 1);

error = 0;

%%
[val, grad_y] = obj.g(y, t);
[val_AD, grad_y_AD] = obj_AD.g(y, t);
local_error = norm(val - val_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end
local_error = norm(grad_y - grad_y_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
[val, grad_z] = obj.R(z);
[val_AD, grad_z_AD] = obj_AD.R(z);
local_error = norm(val - val_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end
local_error = norm(grad_z - grad_z_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = obj.g_yy_Apply(v, y, t);
Mv_AD = obj_AD.g_yy_Apply(v, y, t);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = obj.R_zz_Apply(v, z);
Mv_AD = obj_AD.R_zz_Apply(v, z);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
[f, f_y, f_z] = con.f(y, z, t);
[f_AD, f_y_AD, f_z_AD] = con_AD.f(y, z, t);
local_error = norm(f - f_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end
local_error = norm(f_y - f_y_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end
local_error = norm(f_z - f_z_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
[h, hz] = con.h(z);
[h_AD, hz_AD] = con_AD.h(z);
local_error = norm(h - h_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end
local_error = norm(hz - hz_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = con.f_yy_Apply(v, y, z, t, lambda);
Mv_AD = con_AD.f_yy_Apply(v, y, z, t, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = con.f_yz_Apply(v, y, z, t, lambda);
Mv_AD = con_AD.f_yz_Apply(v, y, z, t, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = con.f_zy_Apply(v, y, z, t, lambda);
Mv_AD = con_AD.f_zy_Apply(v, y, z, t, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = con.f_zz_Apply(v, y, z, t, lambda);
Mv_AD = con_AD.f_zz_Apply(v, y, z, t, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = con.h_zz_Apply(v, z, lambda);
Mv_AD = con_AD.h_zz_Apply(v, z, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
if error > 1.e-9
    fprintf(2,'\nautomatic_differentiation/Transient_Thermal failed.\n');
else
    fprintf(1,'\nautomatic_differentiation/Transient_Thermal passed.\n');
end

con_AD.Clear_AD();
