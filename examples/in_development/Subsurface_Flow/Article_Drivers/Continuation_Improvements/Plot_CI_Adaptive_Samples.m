%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Experiment 2b (Driver_CI_Adaptive_Samples), read from its saved
% results.  With Newton polishing that stops on a gradient tolerance, the
% endpoint gradient is small by construction, so the plots compare the cost to
% reach the tolerance and the fraction of samples that reach it.

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-adaptive-samples-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;

cost_label = sprintf('Cost per sample (continuation steps + %.2g \\times Newton steps)', res.settings.newton_cost);
header = sprintf('Mean-to-sample, Newton until ||g|| \\leq %.2g (\\epsilon_{curv} = %g)', ...
                 res.settings.newton_gradient_tol, res.settings.sample_epsilon_curv);

%% 1-2. Cost to reach the tolerance: fixed N_c vs heuristic settings
disp('Plots 1-2: cost to reach the Newton tolerance');
fixed_summaries = cellfun(@(s) s.summary, res.sample_fixed, 'UniformOutput', false);
adapt_summaries = cellfun(@(a) a.summary, res.sample_adapt, 'UniformOutput', false);
summaries = [fixed_summaries{:}, adapt_summaries{:}, res.sample_backtrack.summary];
Plot_CI_Cost_Panels(summaries, header, cost_label);

%% 3. Does the heuristic assign more continuation steps to harder samples?
% Difficulty proxy: Newton steps needed after a single fixed continuation step.
disp('Plot 3: assigned continuation steps vs difficulty');
j1 = find(res.settings.sample_fixed_sweep == 1, 1);
if ~isempty(j1)
    difficulty = res.sample_fixed{j1}.newton_steps(:);
    eps_list = res.settings.sample_epsilon_nl;
    num_series = length(eps_list);
    figure;
    hold on;
    for kk = 1:num_series
        assigned = res.sample_adapt{kk}.steps(:);
        offset = 0.12 * (kk - (num_series + 1) / 2);      % separate overlapping integer points
        plot(difficulty + offset, assigned, 'o');
        ok = isfinite(difficulty) & isfinite(assigned);
        if nnz(ok) >= 3 && numel(unique(difficulty(ok))) > 1 && numel(unique(assigned(ok))) > 1
            rho = corr(difficulty(ok), assigned(ok), 'Type', 'Spearman');
        else
            rho = NaN;
        end
        fprintf('eps_nl = %g: Spearman(Newton steps at N_c = 1, assigned continuation steps) = %.3f\n', eps_list(kk), rho);
    end
    xlabel('Newton steps needed after a single continuation step');
    ylabel('Continuation steps chosen by the heuristic');
    title('Heuristic step count vs sample difficulty');
    legend(arrayfun(@(e) sprintf('\\epsilon_{nl} = %g', e), eps_list, 'UniformOutput', false), 'Location', 'best');
    drawnow;
end

disp('Adaptive-samples plots done.');
