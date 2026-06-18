%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BF_Sol_Op_Interface < handle

    properties

    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        [u] = State_Solve(this, z)

        [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)

    end

    %% Constructor
    methods

        function this = BF_Sol_Op_Interface()

        end

    end

end
