classdef MD_Opt_Prob_Interface_Py < MD_Opt_Prob_Interface

    properties
        opt_prob_interface_py
    end

    methods

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            z_out = this.opt_prob_interface_py.Apply_Solution_Operator_z_Jacobian_Transpose(u_in, z);
            z_out = double(z_out);
            if size(z_out, 1) == 1
                z_out = z_out';
            end
        end

        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            z_out = this.opt_prob_interface_py.Apply_RS_Hessian(z_in, z);
            z_out = double(z_out);
            if size(z_out, 1) == 1
                z_out = z_out';
            end
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            grad_u = this.opt_prob_interface_py.Misfit_Gradient(u, z);
            grad_u = double(grad_u);
            if size(grad_u, 1) == 1
                grad_u = grad_u';
            end
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = this.opt_prob_interface_py.Apply_Misfit_Hessian(u_in, u, z);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function [u_out] = State_Solve(this, z)
            u_out = this.opt_prob_interface_py.State_Solve(z);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function [u_out] = Apply_Solution_Operator_z_Jacobian(this, z_in, z)
            u_out = this.opt_prob_interface_py.Apply_Solution_Operator_z_Jacobian(z_in, z);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function this = MD_Opt_Prob_Interface_Py(opt_prob_interface_py)
            this.opt_prob_interface_py = opt_prob_interface_py;
        end

    end

end
