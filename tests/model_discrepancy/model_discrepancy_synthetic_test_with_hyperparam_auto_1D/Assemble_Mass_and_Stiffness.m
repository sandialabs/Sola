function [M,S] = Assemble_Mass_and_Stiffness(m)
x = linspace(0, 1, m)';
h = x(2) - x(1);
M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
M(1, 1) = .5 * M(1, 1);
M(end, end) = .5 * M(end, end);
M = (1 / 6) * h * M;
S = diag(2 * ones(1, m)) + (-1) * diag(ones(1, m - 1), 1) + (-1) * diag(ones(1, m - 1), -1);
S(1, 1) = .5 * S(1, 1);
S(end, end) = .5 * S(end, end);
S = (1 / h) * S;
end