%% Clear workspace and add path.
clear;
close all;
clc;
addpath('../../src/model_reduction/operators/');
rng(1120);

%% Set test parameters and generate data.

r = 10;                 % Reduced state dimension.
n_y = 2 * r;            % Full state dimension.
n_q = floor(r / 2);     % Input dimension.
n_t = 4 * n_y;          % Number of snapshots.

Y = randn(n_y, n_t);    % State snapshots.
Q = randn(n_q, n_t);    % Inputs.
Vr = randn(n_y, r);     % Basis matrix.

y = Y(:, 1);            % Single snapshot.
q = Q(:, 1);            % Single input vector.
v = randn(n_y, 1);      % Single search direction.

%% Test Constant_Operator.

op = Constant_Operator();
assert(size(op.entries, 1) == 0);
assert(size(op.entries, 2) == 0);

C = randn(n_y, 1);
op.Set_Entries(C);
assert(all(op.entries == C));

Cy = op.Apply(y);
assert(size(Cy, 1) == n_y);
assert(size(Cy, 2) == 1);
assert(all(Cy == C, 'all'));

CY = op.Apply(Y);
assert(size(CY, 1) == n_y);
assert(size(CY, 2) == n_t);
assert(all(CY(:, 1) == C));

assert(Constant_Operator.Column_Dimension() == 1);

op.Finite_Difference_Check(n_t);
D = Constant_Operator.Datablock(Y);
assert(size(D, 1) == 1);
assert(size(D, 2) == n_t);
assert(all(D == 1));

%% Test Linear_Operator.

A = randn(n_y, n_y);
op = Linear_Operator(A);
assert(all(op.entries == A, 'all'));

AY = op.Apply(Y);
assert(size(AY, 1) == n_y);
assert(size(AY, 2) == n_t);
assert(all(AY == (A * Y), 'all'));

AJac = op.Jacobian_y(Y, Q);
assert(size(AJac, 1) == n_y);
assert(size(AJac, 2) == n_y);
assert(all(AJac == A, 'all'));

assert(Linear_Operator.Column_Dimension(r) == r);

op.Finite_Difference_Check(n_t);
D = Linear_Operator.Datablock(Y);
assert(size(D, 1) == n_y);
assert(size(D, 2) == n_t);
assert(all(D == Y, 'all'));

%% Test Quadratic_Operator.

H = randn(n_y, n_y^2);
op = Quadratic_Operator(H);
assert(size(op.entries, 1) == n_y);
assert(size(op.entries, 2) == n_y * (n_y + 1) / 2);

y = Y(:, 1);
Hy = op.Apply(y);
assert(size(Hy, 1) == n_y);
assert(size(Hy, 2) == 1);
assert(allclose(Hy, H * kron(y, y)));

HY = op.Apply(Y);
assert(size(AY, 1) == n_y);
assert(size(AY, 2) == n_t);
HYtrue = zeros(n_y, n_t);
for j = 1:n_t
    HYtrue(:, j) = H * kron(Y(:, j), Y(:, j));
end
assert(allclose(HY, HYtrue));

I = eye(n_y);
for j = 1:n_t
    y = Y(:, j);
    HJac = op.Jacobian_y(y);
    assert(size(HJac, 1) == n_y);
    assert(size(HJac, 2) == n_y);
    HJactrue = H * (kron(I, y) + kron(y, I));
    assert(allclose(HJac, HJactrue));
end

HHess_v = op.Hessian_yy_Apply(v, y, q, y);
assert(size(HHess_v, 1) == n_y);
assert(size(HHess_v, 2) == size(v, 2));
HHess_V = op.Hessian_yy_Apply(Vr, y, q, y);
assert(size(HHess_V, 1) == n_y);
assert(size(HHess_V, 2) == size(Vr, 2));

assert(Quadratic_Operator.Column_Dimension(r) == r * (r + 1) / 2);

op.Finite_Difference_Check(n_t);
D = Quadratic_Operator.Datablock(Y);
assert(size(D, 1) == n_y * (n_y + 1) / 2);
assert(size(D, 2) == n_t);
assert(allclose(op.Apply(Y), HYtrue));

%% Test Input_Operator

B = randn(n_y, n_q);
op = Input_Operator(B);
assert(all(op.entries == B, 'all'));

BQ = op.Apply(Y, Q);
assert(size(BQ, 1) == n_y);
assert(size(BQ, 2) == n_t);
assert(all(BQ == (B * Q), 'all'));

BJac = op.Jacobian_q(Y, Q);
assert(size(BJac, 1) == n_y);
assert(size(BJac, 2) == n_q);
assert(all(BJac == B, 'all'));

assert(Input_Operator.Column_Dimension(r, n_q) == n_q);

op.Finite_Difference_Check(n_t);
D = Input_Operator.Datablock(Y, Q);
assert(size(D, 1) == n_q);
assert(size(D, 2) == n_t);
assert(all(D == Q, 'all'));

%% Helper functions

function [result] = allclose(A1, A2)
    % Return true if ||A1 - A2|| < the tolerance ``TOL``.
    result = (norm(A1 - A2) < 1e-12);
end
