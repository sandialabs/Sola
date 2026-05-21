%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Linear_OED_A_Opt < handle

    % We assume a linear Bayesian inverse problem with a mean zero Gaussian noise
    % model with covariance sigma^2*I, a Gaussian prior, and a linear observation operator

    properties
        likelihood
        inf_dim_prior
        con

        u_dim
        z_dim
        d_dim
        sigma_sq
        trace_samples
        num_trace_samples
        reguarlization_coeff

        forward_operator_sing_vecs_input
        forward_operator_sing_vecs_output
        forward_operator_sing_vals
        forward_operator_rank

        opt_tol
        fun_tol
        iteration_limit
        step_tol
        verbose
    end

    methods (Access = public)

        function this = Linear_OED_A_Opt(likelihood, inf_dim_prior, con, num_trace_samples, reguarlization_coeff)
            this.likelihood = likelihood;
            this.inf_dim_prior = inf_dim_prior;
            this.con = con;
            z = inf_dim_prior.Get_Prior_Mean();
            this.z_dim = length(z);
            u = con.c_z_Apply(z);
            this.u_dim = length(u);
            d = likelihood.Observation_Operator_Apply(u);
            this.d_dim = length(d);
            tmp = this.likelihood.Noise_Precision_Apply(ones(this.d_dim, 1));
            if max(tmp) - min(tmp) > 0
                disp('Error: the noise covariance must be a scalar multiple of the identity');
            end
            this.sigma_sq = 1 / mean(tmp);

            mass_mat_sqrt = Mass_Matrix_Sqrt(inf_dim_prior);
            omega = randn(this.z_dim, num_trace_samples);
            this.trace_samples = mass_mat_sqrt.Matrix_Sqrt_Apply(omega);
            this.num_trace_samples = num_trace_samples;
            this.reguarlization_coeff = reguarlization_coeff;

            this.opt_tol = 10^-8;
            this.fun_tol = 10^-6;
            this.iteration_limit = 10^3;
            this.step_tol = 10^-6;
            this.verbose = true;
        end

        function [sing_vals] = Compute_Forward_Operator_GSVD(this, num_sing_vals, oversampling, num_subspace_iters)
            vec_in = zeros(this.z_dim, 1);
            vec_out = zeros(this.d_dim, 1);
            gsvd = Forward_Operator_GSVD(vec_in, vec_out, this.inf_dim_prior, this.likelihood, this.con);
            [sing_vecs_input, sing_vecs_output, sing_vals] = gsvd.Compute_GSVD(num_sing_vals, oversampling, num_subspace_iters);
            this.forward_operator_sing_vecs_input = sing_vecs_input;
            this.forward_operator_sing_vecs_output = sing_vecs_output;
            this.forward_operator_sing_vals = sing_vals;
            this.forward_operator_rank = length(sing_vals);
        end

        function [evecs, evals] = Compute_Misfit_Hessian_GEVP(this, w)
            gevp = Misfit_Hessian_GEVP(this.forward_operator_sing_vecs_input, this.forward_operator_sing_vecs_output, this.forward_operator_sing_vals, w, this.inf_dim_prior, this.sigma_sq);
            rank = min(length(find(w > 0)), length(this.forward_operator_sing_vals));
            [evecs, evals] = gevp.Compute_GEVP(rank, 0);
        end

        function [Hinv_v] = Compute_Inverse_Hessian_Matvec(this, v, misfit_evecs, misfit_evals)
            tmp1 = this.inf_dim_prior.Laplacian_Like_Transpose_Inverse_Apply(v);
            tmp2 = this.inf_dim_prior.Mass_Matrix_Apply(tmp1) - misfit_evecs * diag(misfit_evals ./ (1 + misfit_evals)) * (misfit_evecs' * tmp1);
            Hinv_v = this.inf_dim_prior.Laplacian_Like_Inverse_Apply(tmp2);
        end

        function [val, grad] = OED_Objective(this, w)
            [val1, grad1] = this.Posterior_Trace_Objective(w);
            [val2, grad2] = this.R(w);
            val = val1 + this.reguarlization_coeff * val2;
            grad = grad1 + this.reguarlization_coeff * grad2;
        end

        function [val, grad] = Posterior_Trace_Objective(this, w)
            [misfit_evecs, misfit_evals] = this.Compute_Misfit_Hessian_GEVP(w);
            Hinv_samples = this.Compute_Inverse_Hessian_Matvec(this.trace_samples, misfit_evecs, misfit_evals);

            tmp1 = this.inf_dim_prior.Laplacian_Like_Apply(Hinv_samples);
            tmp2 = this.inf_dim_prior.Mass_Matrix_Inverse_Apply(tmp1);
            tmp3 = this.inf_dim_prior.Laplacian_Like_Transpose_Apply(tmp2);
            tmp4 = this.forward_operator_sing_vecs_input' * tmp3;
            tmp5 = diag(this.forward_operator_sing_vals) * tmp4;
            tmp6 = this.forward_operator_sing_vecs_output * tmp5;
            Ftilde_E_Hinv_samples = (1 / sqrt(this.sigma_sq)) * tmp6;

            val = (1 / this.num_trace_samples) * trace(this.trace_samples' * Hinv_samples);
            grad = -(1 / this.num_trace_samples) * sum(Ftilde_E_Hinv_samples.^2, 2);
        end

        function [val, grad] = R(this, w)
            val = sum(w);
            grad = 1 + 0 * w;
        end

        function [w] = Optimize_Design(this)
            w0 = rand(this.d_dim, 1);
            lb = zeros(this.d_dim, 1);
            ub = ones(this.d_dim, 1);
            verb = 'iter-detailed';
            if this.verbose == false
                verb = 'none';
            end
            options = optimoptions(@fmincon, ...
                                   'Display', verb, ...
                                   'Algorithm', 'interior-point', ...
                                   'SpecifyObjectiveGradient', true, ...
                                   'OptimalityTolerance', this.opt_tol, ...
                                   'FunctionTolerance', this.fun_tol, ...
                                   'MaxIterations', this.iteration_limit, ...
                                   'StepTolerance', this.step_tol);
            w = fmincon(@(w)this.OED_Objective(w), w0, [], [], [], [], lb, ub, [], options);
        end

    end

end
