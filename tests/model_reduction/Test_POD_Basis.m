%% Clear workspace and add path.
clear;
close all;
clc;
addpath('../../src/model_reduction/');
rng(1117);

%% Get test data.
n_y = 1000;
n_t = 200;
r = 10;
Y = randn(n_y, n_t);

%% Test POD basis without mean shift.

basis = POD_Basis(Y);

assert(size(basis.ybar, 1) == n_y);
assert(size(basis.ybar, 2) == 1);
assert(all(basis.ybar == 0));

assert(size(basis.V, 1) == n_y);
assert(size(basis.V, 2) == n_t);

basis.r = r;
assert(size(basis.V, 1) == n_y);
assert(size(basis.V, 2) == r);

Yhat = basis.Compress(Y);
assert(size(Yhat, 1) == r);
assert(size(Yhat, 2) == n_t);

Yproj = basis.Decompress(Yhat);
assert(size(Yproj, 1) == n_y);
assert(size(Yproj, 2) == n_t);

%% Test POD basis with mean shift.

basis = POD_Basis(Y, true);

assert(size(basis.V, 1) == n_y);
assert(size(basis.V, 2) == n_t);
assert(size(basis.ybar, 1) == n_y);
assert(size(basis.ybar, 2) == 1);
assert(all(basis.ybar == mean(Y, 2)));

basis.r = r;
assert(size(basis.V, 1) == n_y);
assert(size(basis.V, 2) == r);

Yhat = basis.Compress(Y);
assert(size(Yhat, 1) == r);
assert(size(Yhat, 2) == n_t);

Yproj = basis.Decompress(Yhat);
assert(size(Yproj, 1) == n_y);
assert(size(Yproj, 2) == n_t);

%% Test economization

basis = POD_Basis(Y);
basis.economize = false;
basis.r = r;
basis.r = r + 1;
assert(size(basis.V, 2) == r + 1);

basis.economize = true;
basis.r = r;
assert(basis.maxdim == basis.r);
assert(size(basis.singular_vectors, 2) == basis.r);
