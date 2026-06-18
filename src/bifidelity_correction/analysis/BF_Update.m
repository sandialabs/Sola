%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BF_Update < handle

    properties
        sol_op_interface
        opt_prob_interface
    end

    methods

        function this = BF_Update(sol_op_interface, opt_prob_interface)
            arguments
                sol_op_interface BF_Sol_Op_Interface
                opt_prob_interface MD_Opt_Prob_Interface
            end
            this.sol_op_interface = sol_op_interface;
            this.opt_prob_interface = opt_prob_interface;
        end

        function [z_update] = Update(this,z)
            u_lofi = this.opt_prob_interface.State_Solve(z);
            u_hifi = this.sol_op_interface.State_Solve(z);

            J_grad_u = this.opt_prob_interface.Misfit_Gradient(u_lofi,z);
            B = this.sol_op_interface.Apply_Solution_Operator_z_Jacobian_Transpose(J_grad_u,z) - this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(J_grad_u,z);

            discrep = u_hifi - u_lofi;
            Hess_discrep = this.opt_prob_interface.Apply_Misfit_Hessian(discrep,u_lofi,z);
            B = B + this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(Hess_discrep,z);

            z_update = z - this.Apply_RS_Hessian_Inverse(B,z);
        end

        function [z_out] = Apply_RS_Hessian_Inverse(this, z_in, z)
            z_out = 0 * z_in;
            for k = 1:size(z_in, 2)
                tol = 1.e-7;
                max_iter = length(z);
                [z_out(:, k), flag, relres, iter, resvec] = pcg(@(x)this.opt_prob_interface.Apply_RS_Hessian(x, z), z_in(:, k), tol, max_iter);
                if flag ~= 0
                    disp('CG did not converge');
                end
            end
        end

    end

end
