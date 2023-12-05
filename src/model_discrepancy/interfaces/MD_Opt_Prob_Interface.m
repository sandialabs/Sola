classdef MD_Opt_Prob_Interface < handle

    properties

    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)

        [z_out] = Apply_RS_Hessian(this, z_in, z)

        [grad_u] = Misfit_Gradient(this, u, z)

        [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)

    end

    methods

        % This function must be implemented to enable continuation
        function [u] = State_Solve(this, z)
            u = [];
            disp('State_Solve must be implemented to use continuation algorithm');
        end

        % This function must be implemented to enable continuation
        function [u_out] = Apply_Solution_Operator_z_Jacobian(this, z_in, z)
            u_out = [];
            disp('Apply_Solution_Operator_z_Jacobian must be implemented to use continuation algorithm');
        end

        function this = MD_Opt_Prob_Interface()

        end

    end

end
