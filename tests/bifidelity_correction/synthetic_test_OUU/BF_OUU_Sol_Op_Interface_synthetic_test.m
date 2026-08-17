%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BF_OUU_Sol_Op_Interface_synthetic_test < BF_OUU_Sol_Op_Interface

    properties
        Xi
        alpha
    end

    %% Pure virtual functions for user implementation
    methods (Access = public)

        function [u] = State_Solve_Per_Sample(this, z, s)
            u = this.alpha * this.Xi(1,s) * z.^3 + this.Xi(2,s);
        end

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(this, u_in, z, s)
            z_out = this.alpha * this.Xi(1,s) * 3 * diag(z.^2) * u_in;
        end

    end

    %% Constructor
    methods

        function this = BF_OUU_Sol_Op_Interface_synthetic_test(Xi)
            this.Xi = Xi;
            this.alpha = 1.2;
        end

    end

end
