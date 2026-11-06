%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Diagnostic: is the Hessian used in the Newton polishing exact?
%
% The Newton steps converge linearly (about two digits per step) instead of
% quadratically.  A Newton step with an exact Hessian should converge
% quadratically inside the basin, so a Hessian that is accurate to only ~1% would
% explain the observed rate.  This driver checks
%
%   sen_op.Apply_Hessian(v, beta, theta, t)   versus   the central finite difference
%   (g(beta + e v) - g(beta - e v)) / (2 e)   of sen_op.Gradient,
%
% for random unit directions v, at three points of a few samples' mean-to-sample
% runs (all evaluated at t = 1, as in the Newton phase):
%   1. the mean start beta_bar,
%   2. the end of one continuation step (where Newton starts),
%   3. after several Newton steps (near the solution).
% It also reports the symmetry error of the Hessian, the Newton contraction
% actually observed from the same samples, and a check of the base reduced
% Hessian alone (opt_prob_interface.Apply_RS_Hessian versus the finite difference
% of the reduced gradient with zero discrepancy correction) to localize any error.
%
% Reading the numbers: the finite-difference error has a truncation part (~ e^2)
% and a roundoff part (~ eps / e), so the best error over the e sweep is the
% meaningful one.  Best relative error <~ 1e-6: the Hessian is exact.  ~1e-3 to
% 1e-1: the Hessian is approximate (some term is missing or approximated).

ci_num_data_points = 1;                % same problem as the sample experiments
CI_Setup;
CI_Sample_Reference;
clear ci_num_data_points;

check_samples = [1, 5, 12];            % samples to check
num_dirs = 3;                          % random directions per point
num_newton_steps = 5;                  % Newton steps used to reach the 'near solution' point
h_list = 10.^(-2:-1:-8);               % finite-difference steps, relative to max(1, ||beta||)

ob = sen_op.opt_prob_interface;
rng(7);                                % reproducible directions (the breve samples are already built)
V_dirs = randn(r, num_dirs);
V_dirs = V_dirs ./ vecnorm(V_dirs);

point_names = {'mean start', 'after 1 continuation step', sprintf('after %d Newton steps', num_newton_steps)};
num_pts = length(point_names);
all_rel = nan(num_dirs, length(h_list), num_pts, length(check_samples));
best_rel = nan(num_pts, length(check_samples));
asym = nan(num_pts, length(check_samples));
base_best = nan(1, length(check_samples));
contraction = nan(length(check_samples), num_newton_steps);

fprintf('\n============================================================\n');
fprintf('Hessian check (%d data points): Apply_Hessian vs finite differences of Gradient\n', num_data_points);
fprintf('============================================================\n');
fprintf('Newton convergence parameters in use elsewhere: CG tolerance 1e-6 (relative).\n');

for si = 1:length(check_samples)
    idx = check_samples(si);
    theta = MD_Discrepancy_Mean_To_Sample_Trajectory(1, idx);
    time_index = theta.Get_Number_of_Timesteps();      % t = 1

    % Points along the actual solve: continuation (1 step), then Newton steps.
    qn_prec = MD_Quasi_Newton_Preconditioner(md_hessian_analysis);
    pt = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
    [beta_traj, grad_traj] = pt.Pseudo_Time_Continuation_Forward_Euler(theta, num_newton_steps, 0);
    betas = {beta_bar, beta_traj(:, 2), beta_traj(:, end)};
    gn = vecnorm(grad_traj(:, 2:end));                 % ||g|| after continuation and after each Newton step
    contraction(si, :) = gn(2:end) ./ gn(1:end-1);

    fprintf('\n--- sample %d ---\n', idx);
    fprintf('Observed Newton gradient norms (after continuation, then each Newton step):\n  %s\n', sprintf('%.3e  ', gn));
    fprintf('Observed contraction ratios g(k+1)/g(k):\n  %s\n', sprintf('%.3e  ', contraction(si, :)));
    fprintf('%-28s | ||g||      | best rel. err (per direction)        | at h      | Hessian asymmetry\n', 'point');
    for pp = 1:num_pts
        out = Check_Hessian(sen_op, betas{pp}, theta, time_index, V_dirs, h_list);
        all_rel(:, :, pp, si) = out.rel;
        best_rel(pp, si) = median(out.best);
        asym(pp, si) = out.asym;
        g_here = norm(sen_op.Gradient(betas{pp}, theta, time_index));
        fprintf('%-28s | %.3e | %s | %.0e  | %.2e\n', point_names{pp}, g_here, sprintf('%.2e ', out.best), median(h_list(out.best_k)), out.asym);
    end

    % Error vs h for the point where Newton starts (direction 1), to see the FD error curve.
    out = Check_Hessian(sen_op, betas{2}, theta, time_index, V_dirs(:, 1), h_list);
    fprintf('FD error vs h (direction 1, after 1 continuation step):\n  h   : %s\n  err : %s\n', ...
            sprintf('%.0e  ', h_list), sprintf('%.1e  ', out.rel));

    % Base reduced Hessian alone (no discrepancy corrections): Apply_RS_Hessian vs FD of the reduced gradient.
    z0 = z_opt + md_hessian_analysis.Apply_V(betas{2});
    base_rel = nan(num_dirs, length(h_list));
    for j = 1:num_dirs
        vz = md_hessian_analysis.Apply_V(V_dirs(:, j));
        Hz = ob.Apply_RS_Hessian(vz, z0);
        for k = 1:length(h_list)
            e = h_list(k) * max(1, norm(z0)) / norm(vz);
            Hz_fd = (Reduced_Gradient(ob, z0 + e * vz) - Reduced_Gradient(ob, z0 - e * vz)) / (2 * e);
            base_rel(j, k) = norm(Hz - Hz_fd) / norm(Hz_fd);
        end
    end
    base_best(si) = median(min(base_rel, [], 2));
    fprintf('Base reduced Hessian (Apply_RS_Hessian vs FD of reduced gradient), best rel. err (median over directions): %.2e\n', base_best(si));
end

%% Summary
fprintf('\n============================================================\n');
fprintf('Summary (median over checked samples of the best relative error)\n');
fprintf('============================================================\n');
for pp = 1:num_pts
    fprintf('  Apply_Hessian vs FD, %-28s: %.2e   (asymmetry %.2e)\n', point_names{pp}, median(best_rel(pp, :)), median(asym(pp, :)));
end
fprintf('  Base reduced Hessian vs FD (no corrections)       : %.2e\n', median(base_best));
fprintf('  Observed Newton contraction (median over steps)  : %.2e\n', median(contraction(:), 'omitnan'));

overall = median(best_rel(:), 'omitnan');
fprintf('\nReading:\n');
if overall < 1.e-5
    fprintf('  The Hessian matches finite differences to %.1e: Apply_Hessian is exact to the FD accuracy, so the\n', overall);
    fprintf('  linear rate is NOT caused by an inexact Hessian.  Look at the linear solve (CG tolerance, preconditioner)\n');
    fprintf('  or the accuracy of the state solve / gradient instead.\n');
elseif overall < 1.e-3
    fprintf('  Relative Hessian error %.1e: small but could limit the rate at the gradient levels reached.\n', overall);
else
    fprintf('  Relative Hessian error %.1e: Apply_Hessian is only an approximate Hessian, consistent with the\n', overall);
    fprintf('  observed linear convergence (contraction %.1e per step).  ', median(contraction(:), 'omitnan'));
    if median(base_best) < 1.e-5
        fprintf('The base reduced Hessian is exact, so the\n  missing or approximate term is in the discrepancy corrections of Apply_Hessian.\n');
    else
        fprintf('The base reduced Hessian is also\n  inexact (%.1e), so check Apply_RS_Hessian itself (e.g., a Gauss-Newton approximation).\n', median(base_best));
    end
end
fprintf('(For a Newton step with Hessian error eps_H, the expected linear contraction is about eps_H.)\n');

results = struct('check_samples', check_samples, 'h_list', h_list, 'point_names', {point_names}, ...
                 'all_rel', all_rel, 'best_rel', best_rel, 'asym', asym, 'base_best', base_best, ...
                 'contraction', contraction);
results_file = fullfile(ci_dir, 'ci-hessian-check-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

%% Plot: finite-difference error vs step size
disp('Plot: Hessian relative error vs finite-difference step');
figure;
hold on;
for pp = 1:num_pts
    curve = squeeze(median(median(all_rel(:, :, pp, :), 1, 'omitnan'), 4, 'omitnan'));
    loglog(h_list, curve, 'o-');
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Finite-difference step h (relative)');
ylabel('||H v - FD|| / ||FD||  (median)');
title('Apply\_Hessian vs finite differences of Gradient');
legend(point_names, 'Location', 'best');
drawnow;

%% Local helpers
function out = Check_Hessian(sen_op, beta, theta, time_index, V_dirs, h_list)
    nv = size(V_dirs, 2);
    nh = length(h_list);
    rel = nan(nv, nh);
    scale = max(1, norm(beta));
    for j = 1:nv
        v = V_dirs(:, j);
        Hv = sen_op.Apply_Hessian(v, beta, theta, time_index);
        for k = 1:nh
            e = h_list(k) * scale;
            gp = sen_op.Gradient(beta + e * v, theta, time_index);
            gm = sen_op.Gradient(beta - e * v, theta, time_index);
            Hv_fd = (gp - gm) / (2 * e);
            rel(j, k) = norm(Hv - Hv_fd) / norm(Hv_fd);
        end
    end
    out.rel = rel;
    [out.best, out.best_k] = min(rel, [], 2);
    if nv >= 2
        v = V_dirs(:, 1);
        w = V_dirs(:, 2);
        Hv = sen_op.Apply_Hessian(v, beta, theta, time_index);
        Hw = sen_op.Apply_Hessian(w, beta, theta, time_index);
        out.asym = abs(v' * Hw - w' * Hv) / (abs(v' * Hw) + abs(w' * Hv));
    else
        out.asym = NaN;
    end
end

function g = Reduced_Gradient(ob, z)
    % Reduced gradient of J(S(z), z) with no discrepancy (same expression as the
    % first two lines of sen_op.Gradient with delta = 0).
    u = ob.State_Solve(z);
    [~, grad_u, grad_z] = ob.Objective_Function(u, z);
    g = grad_z + ob.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z);
end
