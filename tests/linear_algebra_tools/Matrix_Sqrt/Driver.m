clear
close all

m = 100;
mat_sqrt_test = Matrix_Sqrt_Test(m);
M = mat_sqrt_test.M;

v = randn(m,1);
Msqrt_v = mat_sqrt_test.Matrix_Sqrt_Apply(v);
M_v = mat_sqrt_test.Matrix_Sqrt_Apply(Msqrt_v);

error = norm(M_v - M*v)/norm(M*v);

if error > 1.e-8
    disp(error)
end