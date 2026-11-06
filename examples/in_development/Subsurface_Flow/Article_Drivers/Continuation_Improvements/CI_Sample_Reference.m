%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference solve for the posterior samples: mean-to-sample with
% reference_steps steps from the fixed-step mean (CI_Setup must have run).
% This is the slowest part of the experiments, so the result is cached in
% continuation-improvements-sample-reference.mat and reused whenever the
% settings that determine it are unchanged (sampling is seeded by OED_Setup,
% so the same settings give the same samples).  Delete the file to force a
% recompute.

ref_file = fullfile(ci_dir, ['continuation-improvements-sample-reference', cache_tag, '.mat']);
ref_key = struct('num_data_points', num_data_points, 'num_samples', num_samples, ...
                 'fixed_mean_steps', fixed_mean_steps, 'reference_steps', reference_steps, 'reference_newton_steps', reference_newton_steps, ...
                 'n', numel(z_opt));

use_cache = false;
if isfile(ref_file)
    cached = load(ref_file);
    use_cache = isfield(cached, 'ref_key') && isequal(cached.ref_key, ref_key);
end

if use_cache
    sample_ref = cached.sample_ref;
    fprintf('Sample reference loaded from cache (%s).\n', ref_file);
else
    fprintf('Sample reference: mean-to-sample, N_c = %d + %d Newton steps (%d samples)...\n', reference_steps, reference_newton_steps, num_samples);
    sample_ref = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', reference_steps, beta_bar, [], [], sample_indices, reference_newton_steps);
    save(ref_file, 'sample_ref', 'ref_key');
end
fprintf('  reference median/max endpoint ||grad|| = %.3e / %.3e\n', ...
        median(sample_ref.grad_norm, 'omitnan'), max(sample_ref.grad_norm, [], 'omitnan'));
