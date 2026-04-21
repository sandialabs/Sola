%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Quasi_Newton_Preconditioner_Bayesian_Inversion < Quasi_Newton_Preconditioner

    properties
        bayes_inv
        gevp
        lis_evecs
        lis_evals
    end

    methods (Access = public)

        % Overloading with an initialization that exploits structure
        function [z_out] = Apply_Initial_Inverse_Hessian_Approximation(this, z_in)
            z_out = this.bayes_inv.prior.Prior_Covariance_Apply(z_in);
            if ~isempty(this.lis_evals)
                tmp = this.lis_evecs' * z_in;
                tmp = tmp .* this.lis_evals ./ (1 + this.lis_evals);
                z_out = z_out - this.lis_evecs * tmp;
            end
        end

        function [u, z, lambda, theta] = Compute_Nominal_Hessian(this, rank, oversampling)
            [this.lis_evecs, this.lis_evals, u, z, lambda, theta] = this.gevp.Compute_Hessian_GEVP(rank, oversampling);
        end

        function this = Quasi_Newton_Preconditioner_Bayesian_Inversion(z_bar, theta_bar, bayes_inv)
            this@Quasi_Newton_Preconditioner();
            this.bayes_inv = bayes_inv;
            this.gevp = Bayesian_Inversion_Hessian_GEVP(z_bar, theta_bar, bayes_inv);
        end

    end
end
