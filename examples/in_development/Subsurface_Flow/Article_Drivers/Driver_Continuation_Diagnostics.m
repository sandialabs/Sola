%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the OED setup and previously generated OED samples
OED_Setup;
oed_results = load('oed-results.mat');
num_data_points = 3;
data_interface.Set_Z_and_D(oed_results.Z_oed(:, 1:num_data_points), oed_results.D_oed(:, 1:num_data_points));

% Posterior sampling
num_samples = 100;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, true);

% Continuation
num_continuation_steps = 7;
plot_ptile = 50;

%% Hybrid-vs-direct continuation sample diagnostics
%
% This uses the same lazy posterior discrepancy realization for the hybrid
% and direct-continuation endpoints of each checked sample.  It answers:
%   * Are direct sample continuations converged?
%   * Is the hybrid tangent step much larger than the nonlinear continuation
%     displacement?
%   * Does increasing sample-continuation steps move the direct samples?

run_sample_map_diagnostics = true;

if run_sample_map_diagnostics
     fprintf('\n============================================================\n');
     fprintf('Hybrid-vs-direct sample map diagnostics\n');
     fprintf('============================================================\n');

     sample_diag_indices = 1:min(5, num_samples);
     sample_diag_steps = unique([num_continuation_steps, 2 * num_continuation_steps, 30]);

     if isempty(md_hessian_analysis.evals)
          r = length(data_interface.z_opt);
     else
          r = length(md_hessian_analysis.evals);
     end

     % Build one mean continuation point and one dense mean Hessian.  The
     % same sensitivity operator is then reused so each sample index has the
     % same lazy breve realization for hybrid and direct continuation.
     sen_cmp = MD_Continuation_Sensitivity_Operators(md_post_sampling, md_hessian_analysis, false);
     qn_mean_cmp = MD_Quasi_Newton_Preconditioner(md_hessian_analysis);
     pt_mean_cmp = Pseudo_Time_Continuation(zeros(r, 1), sen_cmp, qn_mean_cmp);
     theta_mean_cmp = MD_Discrepancy_Parameter_Trajectory(num_continuation_steps, 0);
     beta_mean_traj_cmp = pt_mean_cmp.Pseudo_Time_Continuation_Forward_Euler(theta_mean_cmp);
     beta_bar_cmp = beta_mean_traj_cmp(:, end);
     time_index_mean_cmp = theta_mean_cmp.Get_Number_of_Timesteps();
     [g_mean_cmp, ~] = sen_cmp.Gradient(beta_bar_cmp, theta_mean_cmp, time_index_mean_cmp);

     H_cmp = zeros(r, r);
     I_r = eye(r);
     for kk = 1:r
          H_cmp(:, kk) = sen_cmp.Apply_Hessian(I_r(:, kk), beta_bar_cmp, theta_mean_cmp, time_index_mean_cmp);
     end
     H_cmp = 0.5 * (H_cmp + H_cmp');

     fprintf('Mean point for sample-map diagnostics uses N_c = %d\n', num_continuation_steps);
     fprintf('Mean ||grad_beta J|| = %.3e\n', norm(g_mean_cmp));
     fprintf('Mean Hessian cond    = %.3e\n\n', cond(H_cmp));

     for sample_idx = sample_diag_indices
          B_k = sen_cmp.Apply_B_hybrid(beta_bar_cmp, theta_mean_cmp, time_index_mean_cmp, sample_idx);
          dbeta_hyb = H_cmp \ B_k;
          beta_hyb = beta_bar_cmp - dbeta_hyb;
          z_hyb = data_interface.z_opt + md_hessian_analysis.Apply_V(beta_hyb);

          fprintf('Sample %d\n', sample_idx);
          fprintf('  ||beta_hyb - beta_bar|| = %.3e\n', norm(beta_hyb - beta_bar_cmp));

          beta_full_by_step = zeros(r, length(sample_diag_steps));
          z_full_by_step = zeros(length(data_interface.z_opt), length(sample_diag_steps));
          grad_full_by_step = zeros(length(sample_diag_steps), 1);

          for jj = 1:length(sample_diag_steps)
               N_samp = sample_diag_steps(jj);
               theta_sample_cmp = MD_Discrepancy_Parameter_Trajectory(N_samp, sample_idx);
               qn_sample_cmp = MD_Quasi_Newton_Preconditioner(md_hessian_analysis);
               pt_sample_cmp = Pseudo_Time_Continuation(zeros(r, 1), sen_cmp, qn_sample_cmp);
               beta_sample_traj_cmp = pt_sample_cmp.Pseudo_Time_Continuation_Forward_Euler(theta_sample_cmp);
               beta_full = beta_sample_traj_cmp(:, end);
               z_full = data_interface.z_opt + md_hessian_analysis.Apply_V(beta_full);

               beta_full_by_step(:, jj) = beta_full;
               z_full_by_step(:, jj) = z_full;

               time_index_sample_cmp = theta_sample_cmp.Get_Number_of_Timesteps();
               [g_sample_cmp, ~] = sen_cmp.Gradient(beta_full, theta_sample_cmp, time_index_sample_cmp);
               grad_full_by_step(jj) = norm(g_sample_cmp);

               fprintf('  direct N_c = %d: ||grad|| = %.3e, ||beta_full-beta_bar|| = %.3e, ||z_full-z_hyb|| = %.3e\n', ...
                       N_samp, grad_full_by_step(jj), norm(beta_full - beta_bar_cmp), norm(z_full - z_hyb));
          end

          z_full_ref = z_full_by_step(:, end);
          for jj = 1:length(sample_diag_steps)
               fprintf('    rel ||z_full(N=%d)-z_full(N=%d)|| = %.3e\n', ...
                       sample_diag_steps(jj), sample_diag_steps(end), ...
                       norm(z_full_by_step(:, jj) - z_full_ref) / max(1, norm(z_full_ref)));
          end

          fprintf('  rel ||z_hyb-z_full(N=%d)|| = %.3e\n\n', ...
                  sample_diag_steps(end), norm(z_hyb - z_full_ref) / max(1, norm(z_full_ref)));
     end
end

%% Hybrid diagnostics: dependence on continuation step count
%
% These diagnostics are meant to separate three possible causes of changing
% hybrid uncertainty bands:
%   (1) the mean continuation point is not converged,
%   (2) the final beta-Hessian solve is inaccurate/sensitive to the QN
%       preconditioner,
%   (3) the optimizer map is too nonlinear for endpoint linearization.

run_hybrid_diagnostics = true;

if run_hybrid_diagnostics
     fprintf('\n============================================================\n');
     fprintf('Hybrid continuation diagnostics\n');
     fprintf('============================================================\n');

     diagnostic_steps = unique([1, 2, 3, 5, 10, 20, num_continuation_steps]);
     num_diag_steps = length(diagnostic_steps);

     if isempty(md_hessian_analysis.evals)
          r = length(data_interface.z_opt);
     else
          r = length(md_hessian_analysis.evals);
     end

     grad_norms = zeros(num_diag_steps, 1);
     beta_norms = zeros(num_diag_steps, 1);
     hess_conds = zeros(num_diag_steps, 1);
     hess_sym_errs = zeros(num_diag_steps, 1);
     pcg_direct_errs = zeros(num_diag_steps, 1);
     direct_band_mean_widths = zeros(num_diag_steps, 1);
     direct_band_max_widths = zeros(num_diag_steps, 1);
     z_bar_diag = zeros(length(z_lofi), num_diag_steps);

     for jj = 1:num_diag_steps
          N_diag = diagnostic_steps(jj);
          fprintf('\nDiagnostics for num_continuation_steps = %d\n', N_diag);
          fprintf('------------------------------------------------------------\n');

          hyb_diag = MD_Hybrid_Continuation_Update(md_post_sampling, md_hessian_analysis, N_diag);
          [~, z_bar, beta_bar] = hyb_diag.Posterior_Update_Mean();
          z_bar_diag(:, jj) = z_bar;

          time_index = hyb_diag.theta_traj.Get_Number_of_Timesteps();
          [g_bar, ~] = hyb_diag.sen_op.Gradient(beta_bar, hyb_diag.theta_traj, time_index);
          grad_norms(jj) = norm(g_bar);
          beta_norms(jj) = norm(beta_bar);

          % Build dense beta-Hessian.  In this example r is small, so this is
          % the simplest way to check conditioning and compare direct vs PCG
          % solves used by the hybrid sampler.
          I_r = eye(r);
          H_beta = zeros(r, r);
          for kk = 1:r
               H_beta(:, kk) = hyb_diag.sen_op.Apply_Hessian(I_r(:, kk), beta_bar, hyb_diag.theta_traj, time_index);
          end

          H_beta_sym = 0.5 * (H_beta + H_beta');
          hess_sym_errs(jj) = norm(H_beta - H_beta', 'fro') / max(1, norm(H_beta, 'fro'));
          hess_conds(jj) = cond(H_beta_sym);
          hess_eigs = eig(H_beta_sym);

          % Direct-solve hybrid samples and PCG-vs-direct solve check.
          z_hyb_direct = zeros(length(z_lofi), num_samples);
          num_solve_checks = min(5, num_samples);
          solve_errs = zeros(num_solve_checks, 1);

          for sample_idx = 1:num_samples
               B_k = hyb_diag.sen_op.Apply_B_hybrid(beta_bar, hyb_diag.theta_traj, time_index, sample_idx);
               dbeta_direct = H_beta_sym \ B_k;
               beta_k_direct = beta_bar - dbeta_direct;
               z_hyb_direct(:, sample_idx) = data_interface.z_opt + md_hessian_analysis.Apply_V(beta_k_direct);

               if sample_idx <= num_solve_checks
                    dbeta_pcg = hyb_diag.pt_cont.Apply_Inverse_Hessian(B_k, beta_bar, hyb_diag.theta_traj, time_index);
                    solve_errs(sample_idx) = norm(dbeta_pcg - dbeta_direct) / max(1, norm(dbeta_direct));
               end
          end

          pcg_direct_errs(jj) = max(solve_errs);

          z_hyb_direct_low = prctile(z_hyb_direct, 50 - plot_ptile / 2, 2);
          z_hyb_direct_high = prctile(z_hyb_direct, 50 + plot_ptile / 2, 2);
          direct_band_width = z_hyb_direct_high - z_hyb_direct_low;
          direct_band_mean_widths(jj) = mean(direct_band_width);
          direct_band_max_widths(jj) = max(direct_band_width);

          fprintf('||grad_beta J(beta_bar, theta_bar)|| = %.3e\n', grad_norms(jj));
          fprintf('||beta_bar||                         = %.3e\n', beta_norms(jj));
          fprintf('Hessian symmetry rel err             = %.3e\n', hess_sym_errs(jj));
          fprintf('eig(H_beta_sym) min/max              = %.3e / %.3e\n', min(hess_eigs), max(hess_eigs));
          fprintf('cond(H_beta_sym)                     = %.3e\n', hess_conds(jj));
          fprintf('max PCG-vs-direct solve rel err      = %.3e\n', pcg_direct_errs(jj));
          fprintf('direct hybrid IQR mean/max width     = %.3e / %.3e\n', ...
                  direct_band_mean_widths(jj), direct_band_max_widths(jj));
     end

     z_ref = z_bar_diag(:, end);
     z_bar_ref_dists = zeros(num_diag_steps, 1);
     for jj = 1:num_diag_steps
          z_bar_ref_dists(jj) = norm(z_bar_diag(:, jj) - z_ref) / max(1, norm(z_ref));
     end

     fprintf('\nSummary relative to N = %d mean continuation point:\n', diagnostic_steps(end));
     fprintf('N\t||grad||\t rel ||z_bar-z_ref||\t cond(H)\t PCG/direct\t mean IQR width\n');
     for jj = 1:num_diag_steps
          fprintf('%d\t%.3e\t %.3e\t\t %.3e\t %.3e\t %.3e\n', ...
                  diagnostic_steps(jj), grad_norms(jj), z_bar_ref_dists(jj), ...
                  hess_conds(jj), pcg_direct_errs(jj), direct_band_mean_widths(jj));
     end

     figure;
     tiledlayout(2, 2);

     nexttile;
     semilogy(diagnostic_steps, grad_norms, 'o-');
     xlabel('$N_c$', 'Interpreter', 'latex');
     ylabel('$\|\nabla_\beta J\|$', 'Interpreter', 'latex');
     title('Mean optimality residual', 'Interpreter', 'latex');

     nexttile;
     semilogy(diagnostic_steps, z_bar_ref_dists, 'o-');
     xlabel('$N_c$', 'Interpreter', 'latex');
     ylabel('Relative distance', 'Interpreter', 'latex');
     title('$z_{\rm bar}$ distance to finest diagnostic mean', 'Interpreter', 'latex');

     nexttile;
     semilogy(diagnostic_steps, hess_conds, 'o-');
     xlabel('$N_c$', 'Interpreter', 'latex');
     ylabel('$\mathrm{cond}(H_\beta)$', 'Interpreter', 'latex');
     title('Hybrid Hessian conditioning', 'Interpreter', 'latex');

     nexttile;
     plot(diagnostic_steps, direct_band_mean_widths, 'o-');
     hold on;
     plot(diagnostic_steps, direct_band_max_widths, 's-');
     xlabel('$N_c$', 'Interpreter', 'latex');
     ylabel('IQR width', 'Interpreter', 'latex');
     legend('mean width', 'max width', 'Location', 'best');
     title('Direct-solve hybrid band width', 'Interpreter', 'latex');
     hold off;
end
