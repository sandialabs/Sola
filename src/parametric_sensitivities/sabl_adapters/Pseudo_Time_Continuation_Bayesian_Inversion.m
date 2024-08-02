classdef Pseudo_Time_Continuation_Bayesian_Inversion < Pseudo_Time_Continuation

    properties
        bayes_inv
        gevp
        lis_evecs
        lis_evals
    end

    methods (Access = public)

        % Overloading with an initialization that exploits structure
        function [z_out] = Apply_Initial_BFGS_Inverse_Hessian(this, z_in)
            z_out = this.bayes_inv.prior.Prior_Covariance_Apply(z_in);
            if ~isempty(this.lis_evals)
                tmp = this.lis_evecs' * z_in;
                tmp = tmp .* this.lis_evals ./ (1 + this.lis_evals);
                z_out = z_out - this.lis_evecs * tmp;
            end
        end

        function [] = Compute_Nominal_Hessian(this, rank, oversampling)
            [this.lis_evecs, this.lis_evals, u, z, lambda, theta] = this.gevp.Compute_Hessian_GEVP(rank, oversampling);
            this.sen_op.current_u = u;
            this.sen_op.current_z = z;
            this.sen_op.current_lambda = lambda;
            this.sen_op.current_theta = theta;
        end

        function this = Pseudo_Time_Continuation_Bayesian_Inversion(z_bar, theta_bar, sen_op, bayes_inv)
            this@Pseudo_Time_Continuation(z_bar, theta_bar, sen_op);
            this.bayes_inv = bayes_inv;
            this.gevp = Bayesian_Inversion_Hessian_GEVP(z_bar, theta_bar, bayes_inv);
        end

    end
end
