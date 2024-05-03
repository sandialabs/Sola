classdef Pseudo_Time_Continuation_BFGS_Example_1 < Pseudo_Time_Continuation_BFGS

    properties
        psen_op_nom
    end

    methods (Access = public)

        function [z_out] = Apply_Nominal_Inv_Hessian(this, z_in)

            tol = 1.e-7;
            max_iter = length(z_in);
            [z_out, flag, relres, iter, resvec] = pcg(@(x)this.psen_op_nom.Apply_RS_Hessian(x), z_in, tol, max_iter);
            if flag ~= 0
                disp('CG did not converge');
            end

        end

    end

    methods

        function this = Pseudo_Time_Continuation_BFGS_Example_1(obj, pcon, z_nom, theta_nom)
            this@Pseudo_Time_Continuation_BFGS(obj, pcon, z_nom, theta_nom);

            this.psen_op_nom = Parameteric_Sensitivity_Operators(obj, pcon);
            this.psen_op_nom.Solve_Forward_and_Adjoint_Problems(z_nom, theta_nom);
        end

    end
end
