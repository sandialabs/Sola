%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Quasi_Newton_Preconditioner < Quasi_Newton_Preconditioner

    properties
        hessian_analysis
    end

    methods (Access = public)

        function [beta_out] = Apply_Initial_Inverse_Hessian_Approximation(this, beta_in)
            if isempty(this.hessian_analysis.evals)
                beta_out = beta_in;
            else
                beta_out = beta_in ./ this.hessian_analysis.evals;
            end
        end

        function this = MD_Quasi_Newton_Preconditioner(hessian_analysis)
            arguments
                hessian_analysis MD_Hessian_Analysis
            end
            this@Quasi_Newton_Preconditioner();
            this.hessian_analysis = hessian_analysis;
        end

    end
end
