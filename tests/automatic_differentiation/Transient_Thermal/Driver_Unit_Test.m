clear;
close all;
run('../../../src/Set_Paths');
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
[val, grad_y] = obj.Time_Instance_Objective(y, t);
[val_AD, grad_y_AD] = obj_AD.Time_Instance_Objective(y, t);
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
[val, grad_z] = obj.Regularization_Objective(z);
[val_AD, grad_z_AD] = obj_AD.Regularization_Objective(z);
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
Mv = obj.Time_Instance_Objective_yy_Apply(v, y, t);
Mv_AD = obj_AD.Time_Instance_Objective_yy_Apply(v, y, t);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = obj.Regularization_Objective_zz_Apply(v, z);
Mv_AD = obj_AD.Regularization_Objective_zz_Apply(v, z);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
[f, f_y, f_z] = con.Time_Instance_RHS(y, z, t);
[f_AD, f_y_AD, f_z_AD] = con_AD.Time_Instance_RHS(y, z, t);
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
[h, hz] = con.Initial_Condition(z);
[h_AD, hz_AD] = con_AD.Initial_Condition(z);
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
Mv = con.Time_Instance_RHS_yy_Apply(v, y, z, t, lambda);
Mv_AD = con_AD.Time_Instance_RHS_yy_Apply(v, y, z, t, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = con.Time_Instance_RHS_yz_Apply(v, y, z, t, lambda);
Mv_AD = con_AD.Time_Instance_RHS_yz_Apply(v, y, z, t, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = con.Time_Instance_RHS_zy_Apply(v, y, z, t, lambda);
Mv_AD = con_AD.Time_Instance_RHS_zy_Apply(v, y, z, t, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = con.Time_Instance_RHS_zz_Apply(v, y, z, t, lambda);
Mv_AD = con_AD.Time_Instance_RHS_zz_Apply(v, y, z, t, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
Mv = con.Initial_Condition_zz_Apply(v, z, lambda);
Mv_AD = con_AD.Initial_Condition_zz_Apply(v, z, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

%%
if error > 1.e-9
    disp('Error in automatic differentiation Transient Thermal example');
end

evalc('rmdir("AdiGator_Files/", "s")');
