%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BF_Sol_Op_Interface_synthetic_test < BF_Sol_Op_Interface

    properties

    end

    %% Pure virtual functions for user implementation
    methods (Access = public)

        function [u] = State_Solve(this, z)
            u = 1.2 * z.^3;
        end

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            z_out = 1.2 * 3 * diag(z.^2) * u_in;
        end

    end

    %% Constructor
    methods

        function this = BF_Sol_Op_Interface_synthetic_test()

        end

    end

end
