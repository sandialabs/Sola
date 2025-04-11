classdef MD_OED_Interface < handle

    properties
        data_interface
    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        % Compute L*v, where Sigma = L*L^T
        [L_v] = Apply_Design_Cov_Factor(this, v)

        % Compute Sigma^{-1}*v
        [Sigmainv_v] = Apply_Design_Cov_Inverse(this, v)

    end

    %% Constructor and helper functions
    methods

        function this = MD_OED_Interface(data_interface)
            arguments
                data_interface {MD_Data_Interface}
            end
            this.data_interface = data_interface;
        end

        function [Z] = Compute_Prior_Samples(this, num_samples)
            Omega = randn(length(this.data_interface.z_opt), num_samples);
            Z = this.data_interface.z_opt + this.Apply_Design_Cov_Factor(Omega);
        end

    end

end
