%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_OUU_u_Prior_Interface < MD_u_Prior_Interface

    properties
        us_prior_interface
        data_interface
        ensemble_weighting
    end

    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
            u_out = 0 * u_in;
            n = size(u_in, 2);
            for k = 1:n
                u_in_k = this.data_interface.Reshape_State_to_Mat(u_in(:, k));
                tmp = this.ensemble_weighting.W_s * this.us_prior_interface.Apply_M_u(u_in_k)';
                u_out(:, k) = this.data_interface.Reshape_State_to_Vec(tmp');
            end
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            u_out = 0 * u_in;
            n = size(u_in, 2);
            for k = 1:n
                u_in_k = this.data_interface.Reshape_State_to_Mat(u_in(:, k));
                tmp = this.ensemble_weighting.W_s_inv * this.us_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u_in_k, scalar)';
                u_out(:, k) = this.data_interface.Reshape_State_to_Vec(tmp');
            end
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = 0 * u_in;
            n = size(u_in, 2);
            for k = 1:n
                u_in_k = this.data_interface.Reshape_State_to_Mat(u_in(:, k));
                tmp = this.ensemble_weighting.W_s_inv * this.us_prior_interface.Apply_W_u_Inverse(u_in_k)';
                u_out(:, k) = this.data_interface.Reshape_State_to_Vec(tmp');
            end
        end

        function [u_out] = Sample_with_Covariance_W_u_Inverse(this, num_samples)
            u_out = zeros(this.data_interface.n_u * this.data_interface.n_r, num_samples);
            for k = 1:num_samples
                u_samps = this.us_prior_interface.Sample_with_Covariance_W_u_Inverse(this.data_interface.n_r);
                u_out(:, k) = this.data_interface.Reshape_State_to_Vec((this.ensemble_weighting.R_inv * u_samps')');
            end
        end

        function [u_out] = Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            u_out = zeros(this.data_interface.n_u * this.data_interface.n_r, num_samples);
            for k = 1:num_samples
                u_samps = this.us_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this.data_interface.n_r, scalar);
                u_out(:, k) = this.data_interface.Reshape_State_to_Vec((this.ensemble_weighting.R_inv * u_samps')');
            end
        end

        function this = MD_OUU_u_Prior_Interface(us_prior_interface, data_interface, ensemble_weighting)
            arguments
                us_prior_interface MD_u_Prior_Interface
                data_interface MD_Data_Interface
                ensemble_weighting MD_OUU_Ensemble_Weighting_Matrix
            end
            this.us_prior_interface = us_prior_interface;
            this.data_interface = data_interface;
            this.ensemble_weighting = ensemble_weighting;
        end

    end

end
