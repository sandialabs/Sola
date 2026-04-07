clear;
close all;

m = 100;
h = 1 / (m - 1);
M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
M(1, 1) = .5 * M(1, 1);
M(end, end) = .5 * M(end, end);
M = (1 / 6) * h * M;
M = sparse(M);
L = ichol(M);

mat_sqrt_test = Sparse_Matrix_Sqrt(M,L);

v = linspace(0,1,m)';
[Msqrt_v,rel_res] = mat_sqrt_test.Matrix_Sqrt_Apply(v);

Linv = L \ eye(m);
S = L * sqrtm(Linv * M * Linv');
v_test = S * v;

error = norm(Msqrt_v - v_test) / norm(v_test);

if error > 1.e-8
    fprintf(2, '\nlinear_algebra_tools/Sparse_Matrix_Sqrt failed.\n');
else
    fprintf(1, '\nlinear_algebra_tools/Sparse_Matrix_Sqrt passed.\n');
end
