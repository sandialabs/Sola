classdef Tracer_LoFi_Constraint < Constraint
    properties
        m
        x
        M_z
        S_z
        M_u
        S_u
    end
    methods (Access = public)

        function this = Tracer_LoFi_Constraint()
            % Call the constructor of the superclass
            this@Constraint();

            % Set Mass and Stiffness matrices
            this.M_z = double(py.fluid_flow_1d_lofi.M_one);
            this.S_z = double(py.fluid_flow_1d_lofi.K_mat_one);
            this.M_u = double(py.fluid_flow_1d_lofi.M);
            this.S_u = double(py.fluid_flow_1d_lofi.K_mat);
            this.x = double(py.fluid_flow_1d_lofi.mesh_coordinates);
            this.m = double(py.fluid_flow_1d_lofi.num_steps) * length(this.x);
        end

        function [u] = State_Solve(this, z)
            % Call the Python function state_solve
            u = py.fluid_flow_1d_lofi.state_solve(z, "vector", return_all = true);
            u = double(u)';
        end

        function [u] = State_Solve_Terminal(this, z)
            % Call the Python function state_solve
            u = py.fluid_flow_1d_lofi.state_solve(z, "vector");
            u = double(u)';
        end

        function [u] = State_Solve_Vertex(this, z)
            % Call the Python function state_solve
            u = py.fluid_flow_1d_lofi.state_solve(z, "vertex");
            u = double(u)';
        end

        function [u_out] = c_u_Transpose_Inverse_Apply(this, u_in, u, z)
            error("Not implemented.");
        end

        function [z_out] = c_z_Transpose_Apply(this, u_in, u, z)
            error("Not implemented.");
        end

        function [u_out] = c_u_Inverse_Apply(this, u_in, u, z)
            error("Not implemented.");
        end

        function [con] = c(this, u, z)
            error("Not implemented.");
        end

        function [z_out] = c_z_Apply(this, z_in, u, z)
            error("Not implemented.");
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            error("Not implemented.");
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            error("Not implemented.");
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            error("Not implemented.");
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            error("Not implemented.");
        end

    end
end
