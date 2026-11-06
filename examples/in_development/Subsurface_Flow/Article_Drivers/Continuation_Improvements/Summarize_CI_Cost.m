%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = Summarize_CI_Cost(run_out, newton_cost, newton_gradient_tol, label)
    % Per-sample cost summary of a Run_Samples output when Newton polishing
    % stops on a gradient tolerance.
    %   cost = continuation steps (+ backtracks) + newton_cost * Newton steps taken
    %   converged = the endpoint gradient reached newton_gradient_tol within the
    %               maximum number of Newton steps (i.e., the continuation put the
    %               sample in Newton's basin).
    s.label = label;
    s.cont_steps = run_out.steps(:) + run_out.backtracks(:);
    s.newton_steps = run_out.newton_steps(:);
    s.cost = s.cont_steps + newton_cost * s.newton_steps;
    s.grad_norm = run_out.grad_norm(:);
    s.converged = isfinite(s.grad_norm) & (s.grad_norm <= newton_gradient_tol * (1 + 1.e-8));
    s.failures = run_out.failures(:);

    ok = isfinite(s.cost);
    s.mean_cont = mean(s.cont_steps, 'omitnan');
    s.mean_newton = mean(s.newton_steps, 'omitnan');
    s.mean_cost = mean(s.cost, 'omitnan');
    s.max_cost = max(s.cost, [], 'omitnan');
    if any(ok)
        s.p90_cost = prctile(s.cost(ok), 90);
    else
        s.p90_cost = NaN;
    end
    s.converged_frac = mean(s.converged);
    s.time_per_sample = run_out.time / numel(s.cost);
end
