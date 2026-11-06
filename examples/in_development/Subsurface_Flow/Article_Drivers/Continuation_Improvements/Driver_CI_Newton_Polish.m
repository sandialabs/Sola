%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 4: continuation steps vs continuation steps + final Newton steps.
%
% Compares spending the budget on continuation steps (predictor + corrector at
% each pseudo-time) against fewer continuation steps followed by M extra Newton
% steps at t = 1 (option num_newton_steps of
% Pseudo_Time_Continuation_Forward_Euler).  Extra Newton steps cost one
% inverse-Hessian solve and one gradient each, roughly half a continuation
% step, so wall time per sample is reported next to accuracy.
%
% Accuracy: relative error against the mean-to-sample reference
% (CI_Sample_Reference) and the endpoint gradient norm ||grad_beta J||.

CI_Setup;
CI_Sample_Reference;

step_counts = [1, 2, 3, 5];         % continuation steps N_c
newton_counts = [0, 1, 2];          % extra Newton steps M at t = 1
methods_to_run = {'fixed_mts'};        % add 'nominal' to repeat for nominal-to-sample

fprintf('\n--- Experiment 4: continuation steps + final Newton steps ---\n');
exp4 = struct();
for mm = 1:length(methods_to_run)
    method = methods_to_run{mm};
    runs = cell(length(step_counts), length(newton_counts));
    fprintf('\nMethod: %s\n', method);
    for ii = 1:length(step_counts)
        for jj = 1:length(newton_counts)
            run_out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, method, step_counts(ii), beta_bar, [], [], sample_indices, newton_counts(jj));
            run_out.ref_err = Rel_Col_Err(run_out.z, sample_ref.z);
            runs{ii, jj} = run_out;
            fprintf('N_c = %2d, M = %d | median/max ref err: %.2e / %.2e | median/max ||grad||: %.2e / %.2e | time/sample: %.2fs\n', ...
                    step_counts(ii), newton_counts(jj), median(run_out.ref_err, 'omitnan'), max(run_out.ref_err), ...
                    median(run_out.grad_norm, 'omitnan'), max(run_out.grad_norm), run_out.time / num_samples);
        end
    end
    exp4.(method) = runs;
end

fprintf('\nHighlighted comparison (%s): N_c = 5 vs N_c = 3 + 2 Newton\n', methods_to_run{1});
a = exp4.(methods_to_run{1}){step_counts == 5, newton_counts == 0};
b = exp4.(methods_to_run{1}){step_counts == 3, newton_counts == 2};
fprintf('  5 steps      : median ||grad|| %.2e, max %.2e, time/sample %.2fs\n', median(a.grad_norm), max(a.grad_norm), a.time / num_samples);
fprintf('  3 steps + 2 N: median ||grad|| %.2e, max %.2e, time/sample %.2fs\n', median(b.grad_norm), max(b.grad_norm), b.time / num_samples);

results = struct();
results.lowfi_grad_norm = lowfi_grad_norm;
results.settings = struct('num_data_points', num_data_points, 'num_samples', num_samples, ...
                          'fixed_mean_steps', fixed_mean_steps, 'reference_steps', reference_steps, ...
                          'step_counts', step_counts, 'newton_counts', newton_counts, ...
                          'methods', {methods_to_run});
results.exp4 = exp4;
results_file = fullfile(ci_dir, 'ci-newton-polish-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Newton_Polish;
