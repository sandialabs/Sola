%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Bayesian_Inversion_Objective < Objective

    properties
        likelihood
        prior
    end

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            [val1, grad_u] = this.likelihood.Misfit(u);
            [val2, grad_z] = this.prior.Regularization(z);
            val = val1 + val2;
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            Mv = this.likelihood.Misfit_HessVec(v);
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            Mv = 0 * u;
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            Mv = 0 * z;
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            Mv = this.prior.Regularization_HessVec(v);
        end

    end

    methods (Access = public)

        function this = Bayesian_Inversion_Objective(likelihood, prior)
            this.likelihood = likelihood;
            this.prior = prior;
        end

    end
end
