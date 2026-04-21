%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

%% Test parameters

plot_approx = false;
plot_results = false;

n_y = 3;
nts = 10.^(1:6);
num_trials = length(nts);
errors_1 = zeros(num_trials, n_y);
errors_2 = zeros(num_trials, n_y);
errors_4 = zeros(num_trials, n_y);
errors_6 = zeros(num_trials, n_y);

%% Experiments.

for i = 1:num_trials
    n_t = nts(i);
    t = linspace(0, 4, n_t);
    t2 = t(2:end);

    % Define truth and its derivative.
    Y = [sin(t.^2); exp(t / 3) - cos(2 * t); t.^2 - t];
    dYdt = [2 * t2 .* cos(t2.^2); (exp(t2 / 3) / 3) + 2 * sin(2 * t2); 2 * t2 - 1];

    % Finite difference estimates.
    operators = {Linear_Operator()};
    rom = OpInf_ROM_Constraint(n_y, 2, t(end), n_t, zeros(n_y, 1), operators);
    dYdt1 = rom.Estimate_State_ddts(Y);
    dYdt2 = rom.Estimate_State_ddts_2ndOrder(Y);
    dYdt4 = rom.Estimate_State_ddts_4thOrder(Y);
    dYdt6 = rom.Estimate_State_ddts_6thOrder(Y);

    % Errors.
    for j = 1:n_y
        errors_1(i, j) = norm(dYdt1(j, :) - dYdt(j, :)) / norm(dYdt(j, :));
        errors_2(i, j) = norm(dYdt2(j, :) - dYdt(j, :)) / norm(dYdt(j, :));
        errors_4(i, j) = norm(dYdt4(j, :) - dYdt(j, :)) / norm(dYdt(j, :));
        errors_6(i, j) = norm(dYdt6(j, :) - dYdt(j, :)) / norm(dYdt(j, :));
    end

    % Plots.
    if plot_approx
        for j = 1:n_y
            figure;
            plot(t2, dYdt(j, :));
            hold on;
            plot(t2, dYdt1(j, :));
            plot(t2, dYdt2(j, :));
            plot(t2, dYdt4(j, :));
            plot(t2, dYdt6(j, :));
            legend({'truth', '1st order', '2nd order', '4th order', '6th order'});
            title('approximations');
        end
    end
end

%% Results.

assert(all(min(errors_1, [], 1) < 1e-2));
assert(all(min(errors_2, [], 1) < 1e-2));
assert(all(min(errors_4, [], 1) < 1e-2));
assert(all(min(errors_6, [], 1) < 1e-2));

if plot_results
    for j = 1:n_y
        figure;
        loglog(nts, errors_1(:, j));
        hold on;
        loglog(nts, errors_2(:, j));
        loglog(nts, errors_4(:, j));
        loglog(nts, errors_6(:, j));
        legend({'1st order', '2nd order', '4th order', '6th order'});
        title('errors');
    end
end
