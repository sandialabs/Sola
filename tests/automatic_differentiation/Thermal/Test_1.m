clear;
close all;
rng(12);

print_output = false;

n_y = 50;
con = Thermal_Constraint(n_y);
con_AD = Thermal_Constraint_AD(n_y);
con_AD.verbose = print_output;
evalc('con_AD.AD_Initialization()');

obj = Thermal_Objective(con);
obj_AD = Thermal_Objective_AD(n_y, n_y, con);
obj_AD.verbose = print_output;
evalc('obj_AD.AD_Initialization()');

u = randn(n_y, 1);
z = randn(n_y, 1);
v = randn(n_y, 1);

error = 0;

[val, grad_u, grad_z] = obj.J(u, z);
[val_AD, grad_u_AD, grad_z_AD] = obj_AD.J(u, z);
local_error = norm(val - val_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end
local_error = norm(grad_u - grad_u_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end
local_error = norm(grad_z - grad_z_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

Mv = obj.J_uu_Apply(v, u, z);
Mv_AD = obj_AD.J_uu_Apply(v, u, z);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

Mv = obj.J_uz_Apply(v, u, z);
Mv_AD = obj_AD.J_uz_Apply(v, u, z);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

Mv = obj.J_zu_Apply(v, u, z);
Mv_AD = obj_AD.J_zu_Apply(v, u, z);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

Mv = obj.J_zz_Apply(v, u, z);
Mv_AD = obj_AD.J_zz_Apply(v, u, z);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

Mv = con.c_z_Apply(v, u, z);
Mv_AD = con_AD.c_z_Apply(v, u, z);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

Mv = con.c_z_Transpose_Apply(v, u, z);
Mv_AD = con_AD.c_z_Transpose_Apply(v, u, z);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

Mv = con.c_u_Transpose_Inverse_Apply(v, u, z);
Mv_AD = con_AD.c_u_Transpose_Inverse_Apply(v, u, z);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

Mv = con.c_u_Inverse_Apply(v, u, z);
Mv_AD = con_AD.c_u_Inverse_Apply(v, u, z);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

u = con.State_Solve(z);
u_AD = con_AD.State_Solve(z);
local_error = norm(u - u_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

lambda = randn(n_y, 1);
Mv = con.c_uu_Apply(v, u, z, lambda);
Mv_AD = con_AD.c_uu_Apply(v, u, z, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

lambda = randn(n_y, 1);
Mv = con.c_uz_Apply(v, u, z, lambda);
Mv_AD = con_AD.c_uz_Apply(v, u, z, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

lambda = randn(n_y, 1);
Mv = con.c_zu_Apply(v, u, z, lambda);
Mv_AD = con_AD.c_zu_Apply(v, u, z, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

lambda = randn(n_y, 1);
Mv = con.c_zz_Apply(v, u, z, lambda);
Mv_AD = con_AD.c_zz_Apply(v, u, z, lambda);
local_error = norm(Mv - Mv_AD);
error = max(error, local_error);
if print_output
    disp(['Error = ', num2str(local_error)]);
end

if error > 1.e-11
    fprintf(2, '\nautomatic_differentiation/Thermal failed.\n');
else
    fprintf(1, '\nautomatic_differentiation/Thermal passed.\n');
end

con_AD.Clear_AD();
