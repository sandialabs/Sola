%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_z_Hyperparameter_Interface_Py < MD_z_Hyperparameter_Interface

    properties
        z_hyperparam_interface_py
    end

    methods

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = this.z_hyperparam_interface_py.Load_Spatial_Node_Data();
            spatial_nodes = double(spatial_nodes);
            if size(spatial_nodes, 2) > 3
                spatial_nodes = spatial_nodes';
            end
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.z_hyperparam_interface_py.Load_Time_Node_Data();
            time_nodes = double(time_nodes);
            if size(time_nodes, 2) > 1
                time_nodes = time_nodes';
            end
        end

        function [u] = State_Solve(this, z)
            m = size(z, 2);
            u = this.z_hyperparam_interface_py.State_Solve(z);
            u = double(u);
            if size(u, 1) == m
                u = u';
            end
        end

        function this = MD_z_Hyperparameter_Interface_Py(z_hyperparam_interface_py, z_type, num_state_solves)
            arguments
                z_hyperparam_interface_py
                z_type string
                num_state_solves = 0
            end
            this@MD_z_Hyperparameter_Interface(z_type, num_state_solves);
            this.z_hyperparam_interface_py = z_hyperparam_interface_py;
        end

    end
end
