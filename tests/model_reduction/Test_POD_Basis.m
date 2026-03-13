%% Clear workspace and add path.
clear;
close all;
addpath('../../src/model_reduction/');
rng(1117);

%% Get test data.
n_y = 1000;
n_t = 200;
r = 10;

Y = randn(n_y, n_t);
w = rand(n_y, 1);
W = diag(w);
Id = eye(r);

%% Test without mean shift or weights.

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

VtrpV = basis.V' * basis.V;
assert(size(VtrpV, 1) == r);
assert(size(VtrpV, 2) == r);
assert(allclose(VtrpV, Id));

%% Test with mean shift, without weights.

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

%% Test without mean shift, with weights.

basis = POD_Basis(Y, false, w);

assert(size(basis.V, 1) == n_y);
assert(size(basis.V, 2) == n_t);

VtrpWV = basis.V' * W * basis.V;
assert(allclose(VtrpWV, eye(n_t)));

basis.r = r;
VtrpWV = basis.V' * W * basis.V;
assert(allclose(VtrpWV, Id));

X = randn(n_y);
W = X' * X + 10 * eye(n_y);
basis = POD_Basis(Y, false, W);

VtrpWV = basis.V' * W * basis.V;
assert(allclose(VtrpWV, eye(n_t)));

basis.r = r;
VtrpWV = basis.V' * W * basis.V;
assert(allclose(VtrpWV, Id));

%% Test without mean shift, with sparse weights.
Wpre = sprand(n_y, n_y, 0.05);
W = ((Wpre + Wpre') / 2) + 10 * speye(n_y);
basis = POD_Basis(Y, false, W, true);

VtrpWV = basis.V' * W * basis.V;
assert(allclose(VtrpWV, eye(n_t)));

basis.r = r;
VtrpWV = basis.V' * W * basis.V;
assert(allclose(VtrpWV, Id));

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

function [result] = allclose(A1, A2)
    % Return true if ||A1 - A2|| / ||A2|| < the tolerance.
    result = (norm(A1 - A2) / norm(A2) < 1e-5);
end
