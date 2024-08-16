classdef MD_Opt_Prob_Interface_Python < MD_Opt_Prob_Interface

    properties
        z_current
        u_current
    end

    methods (Access = public)

        function this = MD_Opt_Prob_Interface_Python(md_interface_data)
            this@MD_Opt_Prob_Interface();
            this.z_current = md_interface_data.z_init;
            this.u_current = md_interface_data.u_init;
        end

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            % Call the Python function eval_c
            z_out = py.fluid_flow_1d_lofi.apply_solution_operator_z_jacobian_transpose(u_in, z);
            z_out = double(z_out);
            if isvector(z_out)
                z_out = z_out';
            end
        end

        function [z_out] = Apply_Solution_Operator_z_Jacobian(this, z_in, z)
            % Call the Python function eval_c
            z_out = py.fluid_flow_1d_lofi.apply_solution_operator_z_jacobian(z_in, z);
            z_out = double(z_out);
            if isvector(z_out)
                z_out = z_out';
            end
        end

        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            z_out = py.fluid_flow_1d_lofi.apply_rs_hessian(z_in, z);
            z_out = double(z_out);
            if isvector(z_out)
                z_out = z_out';
            end
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            grad_u = py.fluid_flow_1d_lofi.misfit_gradient(u, z);
            grad_u = double(grad_u);
            if isvector(grad_u)
                grad_u = grad_u';
            end
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = py.fluid_flow_1d_lofi.apply_misfit_hessian(u_in, u, z);
            u_out = double(u_out);
            if isvector(u_out)
                u_out = u_out';
            end
        end

    end

end
