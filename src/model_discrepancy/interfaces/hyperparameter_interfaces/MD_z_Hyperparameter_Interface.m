classdef MD_z_Hyperparameter_Interface < handle

    properties
        z_type
        num_state_solves
        discrepancy_percent_z_variation

        alpha_z
        beta_z
        beta_t
    end

    %% Virtual functions for user implementation
    methods

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = [];
            disp('Load_Spatial_Node_Data is required for hyperparameter algorithm-based initialization');
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = [];
            disp('Load_Time_Node_Data is required for hyperparameter algorithm-based initialization');
        end

        function [u] = State_Solve(this, z)
            disp('State_Solve is required to estimate alpha_z using low-fidelity solves');
            u = [];
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_z_Hyperparameter_Interface(z_type, num_state_solves)
            arguments
                z_type {string}
                num_state_solves (1, 1) {mustBeNumeric} = 0
            end

            if ~(strcmp(z_type, 'spatial field') || strcmp(z_type, 'transient vector') || strcmp(z_type, 'vector'))
                disp('Error in MD_z_Hyperparameter_Interface: The input z_type should be either "spatial field" "transient vector" or "vector"');
            end

            this.z_type = z_type;
            this.num_state_solves = num_state_solves;
            this.discrepancy_percent_z_variation = 1.0;

            this.alpha_z = 0.0;
            this.beta_z = 0.0;
            this.beta_t = 0.0;
        end

        function [] = Set_alpha_z(this, alpha_z_new)
            this.alpha_z = alpha_z_new;
        end

        function [] = Set_beta_z(this, beta_z_new)
            this.beta_z = beta_z_new;
        end

        function [] = Set_beta_t(this, beta_t_new)
            this.beta_t = beta_t_new;
        end

    end
end
