%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_z_Hyperparameter_Interface_synthetic_test_control_vec < MD_z_Hyperparameter_Interface

    properties
        x
    end

    methods

        function [u] = State_Solve(this, z)
            u = z(1) + z(2) * this.x;
        end

        function this = MD_z_Hyperparameter_Interface_synthetic_test_control_vec(num_state_solves, m)
            this@MD_z_Hyperparameter_Interface('vector', num_state_solves);
            this.x = linspace(0, 1, m)';
        end

    end
end
