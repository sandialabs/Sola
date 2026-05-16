%% Test Rank One Updates Matrix
clear;
close all;
dim = 1000;
A = RankOneUpdatesMatrix(eye(dim));
B = eye(dim);
scale = 5e-2;

samples = [];
for i = 1:20
    z = scale * randn(dim, 1);
    samples(i, :) = z;
    A.Add_Update(z);
    B = B + z * z';
end

% test inverting the matrix
error = zeros(dim, 1);
for i = 1:dim
    x = zeros(dim, 1);
    x(i) = 1;
    error(i) = norm(x - A.Inverse_Apply(B * x));
end

disp(max(error));

% test the determinant update formula
prev_det = log(det(B));

z = scale * randn(dim, 1);


old_gain = log(det(B + z * z')) - prev_det;
fast_gain = log(1 + z' * A.Inverse_Apply(z));

disp(abs(old_gain - fast_gain));