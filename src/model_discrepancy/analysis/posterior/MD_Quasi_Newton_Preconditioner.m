%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Quasi_Newton_Preconditioner < Quasi_Newton_Preconditioner

    properties
        md_hessian_analysis
    end

    methods (Access = public)

        function [z_out] = Apply_Initial_Inverse_Hessian_Approximation(this, z_in)
            if isempty(this.md_hessian_analysis.evals)
                z_out = z_in;
            else
                z_out = z_in ./ this.md_hessian_analysis.evals;
            end
        end

        function this = MD_Quasi_Newton_Preconditioner(md_hessian_analysis)
            this@Quasi_Newton_Preconditioner();
            this.md_hessian_analysis = md_hessian_analysis;
        end

    end
end
