%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BF_OUU_Sol_Op_Interface < handle

    properties

    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        [u] = State_Solve_Per_Sample(this, z, s)

        [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(this, u_in, z, s)

    end

    %% Constructor
    methods

        function this = BF_OUU_Sol_Op_Interface()

        end

    end

end