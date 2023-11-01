classdef Prior_Model_Example_1 < Prior_Model

    properties

    end

    methods (Access = public)

        function [z_out] = Prior_Precision_Apply(this, z_in)
            z_out = z_in;
        end

        function [z_prior_mean] = Get_Prior_Mean(this)
            z_prior_mean = 8 * ones(2, 1);
        end

        function [z_out] = Prior_Covariance_Factor_Apply(this, z_in)
            z_out = z_in;
        end

        function this = Prior_Model_Example_1()

        end

    end

end
