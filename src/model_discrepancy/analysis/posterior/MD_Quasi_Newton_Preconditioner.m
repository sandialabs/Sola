%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Quasi_Newton_Preconditioner < Quasi_Newton_Preconditioner

    properties
        hessian_analysis
    end

    methods (Access = public)

        function [z_out] = Apply_Initial_Inverse_Hessian_Approximation(this, z_in)
            if isempty(this.hessian_analysis.evals)
                z_out = z_in;
            else
                z_out = z_in ./ this.hessian_analysis.evals;
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
