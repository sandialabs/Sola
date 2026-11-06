%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 4b: continuation steps vs continuation steps + final Newton steps
% for the posterior mean (the mean analogue of Driver_CI_Newton_Polish).
% Needs only CI_Setup (no sample solves), so it is fast.  Accuracy is the
% relative error against the fixed-step mean reference and the endpoint
% gradient norm ||grad_beta J||.

CI_Setup;

step_counts = [1, 2, 3, 5, 7, 10];     % continuation steps N_c
newton_counts = [0, 1, 2, 3];             % extra Newton steps M at t = 1

fprintf('\n--- Experiment 4b: posterior mean, continuation steps + final Newton steps ---\n');
runs = cell(length(step_counts), length(newton_counts));
first_run = true;
for ii = 1:length(step_counts)
    for jj = 1:length(newton_counts)
        if first_run
            Run_Mean(sen_op, md_hessian_analysis, z_opt, r, step_counts(ii), [], mean_bt_options, newton_counts(jj));
            first_run = false;
        end
        rr = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, step_counts(ii), [], mean_bt_options, newton_counts(jj));
        rr.ref_err = Rel_Err(rr.z, mean_ref.z);
        runs{ii, jj} = rr;
        fprintf('N_c = %2d, M = %d | ref err %.3e | ||grad|| %.3e | %.2fs\n', ...
                step_counts(ii), newton_counts(jj), rr.ref_err, rr.grad_norm, rr.time);
    end
end

results = struct();
results.settings = struct('num_data_points', num_data_points, 'mean_reference_steps', mean_reference_steps, ...
                          'step_counts', step_counts, 'newton_counts', newton_counts);
results.lowfi_grad_norm = lowfi_grad_norm;
results.mean_ref = mean_ref;  results.runs = runs;
results_file = fullfile(ci_dir, 'ci-mean-newton-polish-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Mean_Newton_Polish;
