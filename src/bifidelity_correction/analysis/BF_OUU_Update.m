%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BF_OUU_Update < handle

    properties
        ouu_sol_op_interface
        ouu_opt_prob_interface
    end

    methods

        function this = BF_OUU_Update(ouu_sol_op_interface, ouu_opt_prob_interface)
            arguments
                ouu_sol_op_interface BF_OUU_Sol_Op_Interface
                ouu_opt_prob_interface MD_OUU_Opt_Prob_Interface
            end
            this.ouu_sol_op_interface = ouu_sol_op_interface;
            this.ouu_opt_prob_interface = ouu_opt_prob_interface;
        end

        function [z_update] = Update(this,z)

            B = 0 * z;
            M = this.ouu_opt_prob_interface.n_r;

            for s = 1:M

                u_lofi = this.ouu_opt_prob_interface.State_Solve_Per_Sample(z, s);
                u_hifi = this.ouu_sol_op_interface.State_Solve_Per_Sample(z, s);

                J_grad_u = this.ouu_opt_prob_interface.ens_weights(s) * this.ouu_opt_prob_interface.Misfit_Gradient_Per_Sample(u_lofi, z, s);
                B = B + this.ouu_sol_op_interface.Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(J_grad_u, z, s) - this.ouu_opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(J_grad_u, z, s);

                discrep = u_hifi - u_lofi;
                Hess_discrep = this.ouu_opt_prob_interface.ens_weights(s) * this.ouu_opt_prob_interface.Apply_Misfit_Hessian_Per_Sample(discrep,u_lofi, z, s);
                B = B + this.ouu_opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(Hess_discrep, z, s);

            end

            z_update = z - this.Apply_RS_Hessian_Inverse(B, z);
        end

        function [z_out] = Apply_RS_Hessian_Inverse(this, z_in, z)
            z_out = 0 * z_in;
            for k = 1:size(z_in, 2)
                tol = 1.e-7;
                max_iter = length(z) + 10;
                [z_out(:, k), flag, relres, iter, resvec] = pcg(@(x)this.ouu_opt_prob_interface.Apply_RS_Hessian(x, z), z_in(:, k), tol, max_iter);
                if flag ~= 0
                    disp('CG did not converge');
                end
            end
        end

    end

end
