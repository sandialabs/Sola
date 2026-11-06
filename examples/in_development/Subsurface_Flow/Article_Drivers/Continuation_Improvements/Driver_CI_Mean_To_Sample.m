%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 1: mean-to-sample vs nominal-to-sample continuation.
%
% Nominal-to-sample starts every sample path at beta = 0 and must traverse
% the whole nonlinear update.  Mean-to-sample starts at the posterior-mean
% update beta_bar (computed once in CI_Setup) and only continues the
% sample-minus-mean discrepancy increment.  Accuracy is measured against the
% mean-to-sample reference (CI_Sample_Reference) and by the endpoint
% optimality residual ||grad_beta J||, which needs no reference.  A third
% variants, mean-to-sample + M extra Newton steps at t = 1 for each M in
% polish_Ms, are run for comparison.

CI_Setup;
CI_Sample_Reference;

step_counts = [1, 2, 3, 4, 5];     % N_c sweep for both sample methods
show_N_c = 2;                          % N_c used by the per-sample plots
polish_Ms = [2, 1];                    % extra Newton steps at t = 1 for mean-to-sample + Newton

fprintf('\n--- Experiment 1: mean-to-sample vs nominal-to-sample ---\n');
num_N = length(step_counts);
exp1 = struct('N_c', num2cell(step_counts), 'nominal', [], 'mts', [], 'mts_newton', []);
for jj = 1:num_N
    N_c = step_counts(jj);
    exp1(jj).nominal = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'nominal', N_c, [], [], [], sample_indices);
    exp1(jj).mts = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', N_c, beta_bar, [], [], sample_indices);
    exp1(jj).mts_newton = cell(1, length(polish_Ms));
    for kk = 1:length(polish_Ms)
        run_out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', N_c, beta_bar, [], [], sample_indices, polish_Ms(kk));
        run_out.ref_err = Rel_Col_Err(run_out.z, sample_ref.z);
        exp1(jj).mts_newton{kk} = run_out;
    end
    exp1(jj).nominal.ref_err = Rel_Col_Err(exp1(jj).nominal.z, sample_ref.z);
    exp1(jj).mts.ref_err = Rel_Col_Err(exp1(jj).mts.z, sample_ref.z);
    fprintf('N_c = %2d | median ref err / ||grad|| / time per sample:\n', N_c);
    fprintf('  nominal : %.2e / %.2e / %.2fs\n', median(exp1(jj).nominal.ref_err, 'omitnan'), ...
            median(exp1(jj).nominal.grad_norm, 'omitnan'), exp1(jj).nominal.time / num_samples);
    fprintf('  MTS     : %.2e / %.2e / %.2fs\n', median(exp1(jj).mts.ref_err, 'omitnan'), ...
            median(exp1(jj).mts.grad_norm, 'omitnan'), exp1(jj).mts.time / num_samples);
    for kk = 1:length(polish_Ms)
        run_out = exp1(jj).mts_newton{kk};
        fprintf('  MTS+%dN  : %.2e / %.2e / %.2fs\n', polish_Ms(kk), median(run_out.ref_err, 'omitnan'), ...
                median(run_out.grad_norm, 'omitnan'), run_out.time / num_samples);
    end
end

results = struct();
results.settings = struct('num_data_points', num_data_points, 'num_samples', num_samples, ...
                          'fixed_mean_steps', fixed_mean_steps, 'reference_steps', reference_steps, ...
                          'step_counts', step_counts, 'show_N_c', show_N_c, 'polish_Ms', polish_Ms);
results.lowfi_grad_norm = lowfi_grad_norm;
results.z_opt = z_opt;  results.z_hifi = z_hifi;
results.mean_start = mean_start;  results.sample_ref = sample_ref;  results.exp1 = exp1;
results_file = fullfile(ci_dir, 'ci-mean-to-sample-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Mean_To_Sample;
