classdef Pseudo_Time_Continuation_BFGS_Bayes_Inv_Prob < Pseudo_Time_Continuation_BFGS

    properties

    end

    methods (Access = public)

        function [z_out] = Apply_Nominal_Inv_Hessian(this, z_in)

        end

    end

    methods

        function this = Pseudo_Time_Continuation_BFGS_Bayes_Inv_Prob(obj, pcon, z_nom, theta_nom)
            this@Pseudo_Time_Continuation_BFGS(obj, pcon, z_nom, theta_nom);
        end

    end
end
