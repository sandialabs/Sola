%% Clear workspace and add path.
clear;
close all;
clc;
addpath('../../src/model_reduction/operators/');
rng(1129);

%% Set test parameters and generate data.

r = 10;                         % Reduced state dimension.
n_y = 2 * r;                    % Full state dimension.
n_q = floor(r / 2);             % Input dimension.
n_t = 4 * n_y;                  % Number of snapshots.
T = 1;                          % Final simulation time.

Y = randn(n_y, n_t);            % State snapshots.
Q = randn(n_q, n_t);            % Inputs.
Vr = randn(n_y, r);             % Basis matrix.

y0 = randn(n_y, 1);             % Initial condition.
y = Y(:, 1);                    % Single snapshot.
q = Q(:, 1);                    % Single input vector.
z = reshape(Q, n_q * n_t, 1);   % Full control vector.

%% Test Constructor.

operators = {Constant_Operator(), Linear_Operator(), ...
             Quadratic_Operator(), Input_Operator()};

con = OpInf_ROM_Constraint(n_y, n_q, T, n_t, y0, operators);

assert(con.n_z == n_q * n_t);
assert(con.n_y == n_y);
assert(con.n_q == n_q);

assert(length(con.ds) == length(operators));
assert(con.ds(1) == 1);
assert(con.ds(2) == n_y);
assert(con.ds(3) == n_y * (n_y + 1) / 2);
assert(con.ds(4) == n_q);
assert(con.d == 1 + n_y + (n_y * (n_y + 1) / 2) + n_q);

%% Test the regularizer.

assert(size(con.regularizer, 1) == con.d);
assert(size(con.regularizer, 2) == con.d);
assert(norm(con.regularizer) == 0);

con.regularizer = 3;
assert(size(con.regularizer, 1) == con.d);
assert(size(con.regularizer, 2) == con.d);
assert(all(diag(con.regularizer) == 3, 'all'));
assert(all(diag(diag(con.regularizer)) == con.regularizer, 'all'));

con.regularizer = [1 2 3 4];
assert(size(con.regularizer, 1) == con.d);
assert(size(con.regularizer, 2) == con.d);
regdiag = diag(con.regularizer);
assert(regdiag(1) == 1);
assert(all(regdiag(2:(n_y + 1)) == 2, 'all'));
assert(all(regdiag((n_y + 2):(n_y * (n_y + 1) / 2 + n_y + 1)) == 3, 'all'));
assert(all(regdiag((n_y * (n_y + 1) / 2 + n_y + 2):end) == 4, 'all'));

%% Test Learn_Operators() (actual Operator Inference part).

con.Learn_Operators(Y, Q);
for i = 1:length(operators)
    entries = con.operators{i}.entries;
    assert(size(entries, 1) > 0);
    assert(min(abs(entries), [], 'all') > 0);
end

%% Finite difference tests.

con.operators = Dummy_Operators(n_y, n_q);
con.verbose = false;

for j = 1:n_t
    t = con.t_mesh(j);
    % disp(['Checks for t = ', num2str(t)]);
    [diffs_y, diffs_z] = con.f_Jacobian_Check(y, z, t);
    checkdiffs(diffs_y);
    checkdiffs(diffs_z);

    [diffs_yy, diffs_yz, diffs_zy, diffs_zz] = con.f_Hessian_Check(y, z, t);
    checkdiffs(diffs_yy);
    checkdiffs(diffs_yz);
    checkdiffs(diffs_zy);
    checkdiffs(diffs_zz);
end

clc;

%% Helper functions

function [result] = allclose(A1, A2)
    % Return true if ||A1 - A2|| < the tolerance ``TOL``.
    result = (norm(A1 - A2) < 1e-12);
end

function checkdiffs(diffs)
    if ~all(isnan(diffs)) && min(diffs > 1e-6)
        error('finite difference check failed');
    end
end

function [operators] = Dummy_Operators(n_y, n_q)
    % Get dummy operators with populated entries.
    op_C = Constant_Operator(randn(n_y, 1));
    op_A = Linear_Operator(randn(n_y, n_y));
    op_H = Quadratic_Operator(randn(n_y, n_y * (n_y + 1) / 2));
    op_B = Input_Operator(randn(n_y, n_q));
    operators = {op_C, op_A, op_H, op_B};
end
