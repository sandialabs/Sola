%% Poisson Test Forward
%
N = 100;
x = linspace(0, 1, N)';
m = 2 * ones(N, 1);

cons = Poisson_Constraint(N, 1);
M = cons.M;
S = cons.S;

true_u = x .* (1 - x);
u = cons.State_Solve(m);

figure;
plot(x, true_u);
hold;
plot(x, u);
legend({'True', 'Solution'})

diff = true_u - u;
disp(diff' * M * diff + diff' * S * diff)
