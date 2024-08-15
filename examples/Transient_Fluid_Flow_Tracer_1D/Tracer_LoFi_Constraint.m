classdef Tracer_LoFi_Constraint < Constraint
    properties
        m
        diff_coeff
        react_coeff
        x
        M
        S
    end
    methods (Access = public)
        function this = Tracer_LoFi_Constraint()
            % Call the constructor of the superclass
            this@Constraint();

            % Set Mass and Stiffness matrices
            this.M = double(py.fluid_flow_1d_lofi.M);
            this.S = double(py.fluid_flow_1d_lofi.K_mat);
            this.m = 31;
        end

        function [u] = State_Solve(this, z)
            % Call the Python function state_solve
            u = py.fluid_flow_1d_lofi.state_solve(z, "vector");
            u = double(u)';
        end

        function [u_out] = c_u_Transpose_Inverse_Apply(this, u_in, u, z)
            % Call the Python function c_u_inv_transpose_apply
            u_out = py.fluid_flow_1d_lofi.c_u_inv_transpose_apply(u_in, z, u);
            u_out = double(u_out)';
        end

        function [z_out] = c_z_Transpose_Apply(this, u_in, u, z)
            % Call the Python function c_z_transpose_apply
            z_out = py.fluid_flow_1d_lofi.c_z_transpose_apply(u_in, z, u);
            z_out = double(z_out)';
        end

        function [u_out] = c_u_Inverse_Apply(this, u_in, u, z)
            % Call the Python function c_u_inv_apply
            u_out = py.fluid_flow_1d_lofi.c_u_inv_apply(u_in, z, u);
            u_out = double(u_out)';
        end

        function [con] = c(this, u, z)
            % Call the Python function eval_c
            con = py.fluid_flow_1d_lofi.eval_c(z, u, return_type="vector");
            con = double(con);
            if isvector(con)
                con = con';
            end
        end

        function [z_out] = c_z_Apply(this, z_in, u, z)
            z_out = py.fluid_flow_1d_lofi.c_z_apply(z_in, z, u);
            z_out = double(z_out)';
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            error("Not implemented.")
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            error("Not implemented.")
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            error("Not implemented.")
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            error("Not implemented.")
        end
    end
end
