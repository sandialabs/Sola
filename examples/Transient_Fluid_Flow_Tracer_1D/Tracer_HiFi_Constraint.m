classdef Tracer_HiFi_Constraint < Constraint
    methods (Access = public)
        function this = Tracer_HiFi_Constraint()
            % Call the constructor of the superclass
            this@Constraint();
        end

        function [u] = State_Solve(this, z)
            % Call the Python function state_solve
            u = py.fluid_flow_1d_hifi_eval.state_solve(z, 'vector');
            u = double(u)'; % Convert Python array to MATLAB array
        end


        function [u_out] = c_u_Transpose_Inverse_Apply(this, u_in, u, z)
            error("Not implemented.")
        end

        function [z_out] = c_z_Transpose_Apply(this, u_in, u, z)
            error("Not implemented.")
        end

        function [u_out] = c_u_Inverse_Apply(this, u_in, u, z)
            error("Not implemented.")
        end

        function [con] = c(this, u, z)
            error("Not implemented.")
        end

        function [z_out] = c_z_Apply(this, u_in, u, z)
            error("Not implemented.")
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
