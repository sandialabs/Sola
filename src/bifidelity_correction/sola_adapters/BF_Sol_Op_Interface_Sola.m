%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BF_Sol_Op_Interface_Sola < BF_Sol_Op_Interface

    properties
        sola_con
        z_current
        u_current
    end

    %% Implementation of base class virtual functions
    methods

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            this.State_Solve(z);
            tmp = this.sola_con.c_u_Transpose_Inverse_Apply(u_in, this.u_current, z);
            z_out = -this.sola_con.c_z_Transpose_Apply(tmp, this.u_current, z);
        end

        function [u] = State_Solve(this, z)
            if norm(z - this.z_current) == 0
                u = this.u_current;
            else
                u = this.sola_con.State_Solve(z);
                this.z_current = z;
                this.u_current = u;
            end
        end

    end

    %% Constructor and helper function
    methods

        function this = BF_Sol_Op_Interface_Sola(sola_con)
            arguments
                sola_con Constraint
            end
            this@BF_Sol_Op_Interface();
            this.sola_con = sola_con;
            this.z_current = inf;
            this.u_current = inf;
        end

    end

end
