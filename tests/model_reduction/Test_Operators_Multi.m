%% Clear workspace and add path.
clear;
close all;
addpath('../../src/model_reduction/operators/');
rng(1120);

enable_fd_test = false;

%% Set test parameters and generate test data.

ns = [3; 6; 2];         % Substate dimensions.
n_y = sum(ns);          % Complete state dimension.
n_q = floor(n_y / 2);   % Input dimension.
n_t = 4 * n_y;          % Number of snapshots.

Y = randn(n_y, n_t);    % State snapshots.
Q = randn(n_q, n_t);    % Inputs.

y = Y(:, 1);            % Single snapshot.
q = Q(:, 1);            % Single input vector.
v = randn(n_y, 1);      % Single search direction.

%% Test Linear_Operator_Multi.

idx = 2;
op = Linear_Operator_Multi(1, idx, ns);
assert(numel(op.entries) == 0);
y_sub = op.Get_Substate(idx, y);
assert(size(y_sub, 1) == ns(idx));
assert(size(y_sub, 2) == 1);
block = op.Datablock(Y);
assert(size(block, 1) == ns(idx));
assert(size(block, 2) == n_t);
assert(op.Column_Dimension() == ns(idx));

A = randn(ns(1), ns(idx));
op.Set_Entries(A);
Ay = op.Apply(y, q);
assert(size(Ay, 1) == ns(1));
assert(size(Ay, 2) == 1);
jac = op.Jacobian_y(y, q);
assert(size(jac, 1) == ns(1));
assert(size(jac, 2) == n_y);
[rows, cols] = find(jac);
minrow = min(rows);
maxrow = max(rows);
mincol = min(cols);
maxcol = max(cols);
jacBlock = jac(minrow:maxrow, mincol:maxcol);
assert(size(jacBlock, 1) == size(A, 1));
assert(size(jacBlock, 2) == size(A, 2));
assert(allclose(jacBlock, A));

assert(op.Jacobian_q(y, q) == 0);
assert(op.Hessian_yy_Apply(y, y, q, y) == 0);
assert(op.Hessian_yq_Apply(q, y, q, y) == 0);
assert(op.Hessian_qy_Apply(y, y, q, y) == 0);
assert(op.Hessian_qq_Apply(q, y, q, y) == 0);

if enable_fd_test
    op.Finite_Difference_Check();
end

%% Test Quadratic_Operator_Multi.

idx1 = 2;
idx2 = 3;
op = Quadratic_Operator_Multi(1, idx1, idx2, ns);
assert(numel(op.entries) == 0);
block = op.Datablock(Y);
assert(size(block, 1) == ns(idx1) * ns(idx2));
assert(size(block, 2) == n_t);
assert(op.Column_Dimension() == ns(idx1) * ns(idx2));

H = randn(ns(1), ns(idx1) * ns(idx2));
op.Set_Entries(H);
Hy = op.Apply(y, q);
assert(size(Hy, 1) == ns(1));
assert(size(Hy, 2) == 1);
jac = op.Jacobian_y(y, q);
assert(size(jac, 1) == ns(1));
assert(size(jac, 2) == n_y);
jacblock1 = op.Get_Substate(idx1, jac')';
assert(size(jacblock1, 1) == ns(1));
assert(size(jacblock1, 2) == ns(idx1));
y2 = op.Get_Substate(idx2, y);
assert(allclose(jacblock1, H * kron(eye(ns(idx1)), y2)));
jacblock2 = op.Get_Substate(idx2, jac')';
y1 = op.Get_Substate(idx1, y);
assert(allclose(jacblock2, H * kron(y1, eye(ns(idx2)))));

if enable_fd_test
    op.Finite_Difference_Check();
end

assert(op.Jacobian_q(y, q) == 0);
assert(op.Hessian_yq_Apply(q, y, q, y) == 0);
assert(op.Hessian_qy_Apply(y, y, q, y) == 0);
assert(op.Hessian_qq_Apply(q, y, q, y) == 0);

%% Test Input_Operator_Multi.

op = Input_Operator_Multi(1, n_q, ns);
assert(numel(op.entries) == 0);
block = op.Datablock(Y, Q);
assert(size(block, 1) == n_q);
assert(size(block, 2) == n_t);
assert(op.Column_Dimension() == n_q);
assert(allclose(block, Q));

B = randn(ns(1), n_q);
op.Set_Entries(B);
Bq = op.Apply(y, q);
assert(size(Bq, 1) == ns(1));
assert(size(Bq, 2) == 1);
jac = op.Jacobian_q(y, q);
assert(size(jac, 1) == ns(1));
assert(size(jac, 2) == n_q);
assert(allclose(jac, B));

assert(op.Jacobian_y(y, q) == 0);
assert(op.Hessian_yy_Apply(y, y, q, y) == 0);
assert(op.Hessian_yq_Apply(q, y, q, y) == 0);
assert(op.Hessian_qy_Apply(y, y, q, y) == 0);
assert(op.Hessian_qq_Apply(q, y, q, y) == 0);

if enable_fd_test
    op.Finite_Difference_Check(30, true);
end

%% Test Input_Squared_Operator_Multi.

op = Input_Squared_Operator_Multi(1, n_q, ns);
assert(numel(op.entries) == 0);
block = op.Datablock(Y, Q);
assert(size(block, 1) == n_q);
assert(size(block, 2) == n_t);
assert(op.Column_Dimension() == n_q);
assert(allclose(block, Q .* Q));

B = randn(ns(1), n_q);
op.Set_Entries(B);
Bq = op.Apply(y, q);
assert(size(Bq, 1) == ns(1));
assert(size(Bq, 2) == 1);
jac = op.Jacobian_q(y, q);
assert(size(jac, 1) == ns(1));
assert(size(jac, 2) == n_q);

assert(op.Jacobian_y(y, q) == 0);
assert(op.Hessian_yy_Apply(y, y, q, y) == 0);
assert(op.Hessian_yq_Apply(q, y, q, y) == 0);
assert(op.Hessian_qy_Apply(y, y, q, y) == 0);
B_qq = op.Hessian_qq_Apply(q, y, q, y);
assert(size(B_qq, 1) == n_q);

if enable_fd_test
    op.Finite_Difference_Check(30, true);
end

%% Helper functions

function [result] = allclose(A1, A2)
    % Return true if ||A1 - A2|| < the tolerance ``TOL``.
    result = (norm(A1 - A2) < 1e-12);
end
