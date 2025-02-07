classdef MD_Hyperparameters_synthetic_test_with_hyperparam_auto < MD_Hyperparameters

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Node_Data(this)
            nodes = this.x;
        end

        function this = MD_Hyperparameters_synthetic_test_with_hyperparam_auto(data_interface,m)
            this@MD_Hyperparameters(data_interface);
            this.x = linspace(0,1,m)';
        end

    end

end