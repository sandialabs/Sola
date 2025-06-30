classdef MD_OUU_u_Prior_Interface < MD_u_Prior_Interface

    properties
        us_prior_interface
        data_interface
        C
        R
        Rinv
    end

    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
            u_out = 0 * u_in;
            n = size(u_in, 2);
            for k = 1:n
                u_in_k = this.data_interface.Reshape_State_to_Mat(u_in(:, k));
                tmp = this.us_prior_interface.Apply_M_u(u_in_k) * this.C;
                u_out(:, k) = this.data_interface.Reshape_State_to_Vec(tmp);
            end
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            u_out = 0 * u_in;
            n = size(u_in, 2);
            for k = 1:n
                u_in_k = this.data_interface.Reshape_State_to_Mat(u_in(:, k));
                tmp = this.us_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u_in_k, scalar) * this.C;
                u_out(:, k) = this.data_interface.Reshape_State_to_Vec(tmp);
            end
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = 0 * u_in;
            n = size(u_in, 2);
            for k = 1:n
                u_in_k = this.data_interface.Reshape_State_to_Mat(u_in(:, k));
                tmp = this.us_prior_interface.Apply_W_u_Inverse(u_in_k) * this.C;
                u_out(:, k) = this.data_interface.Reshape_State_to_Vec(tmp);
            end
        end

        function [u_out] = Sample_with_Covariance_W_u_Inverse(this, num_samples)
            u_samps = this.us_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samples * this.data_interface.n_r);
            u_out = zeros(this.data_interface.n_u * this.data_interface.n_r, num_samples);
            for k = 1:num_samples
                I = ((k - 1) * this.data_interface.n_r + 1):(k * this.data_interface.n_r);
                u_out(:, k) = this.data_interface.Reshape_State_to_Vec(u_samps(:, I) * this.Rinv);
            end
        end

        function [u_out] = Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            u_samps = this.us_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(num_samples * this.data_interface.n_r, scalar);
            u_out = zeros(this.data_interface.n_u * this.data_interface.n_r, num_samples);
            for k = 1:num_samples
                I = ((k - 1) * this.data_interface.n_r + 1):(k * this.data_interface.n_r);
                u_out(:, k) = this.data_interface.Reshape_State_to_Vec(u_samps(:, I) * this.Rinv);
            end
        end

        function this = MD_OUU_u_Prior_Interface(us_prior_interface, data_interface)
            this.us_prior_interface = us_prior_interface;
            this.data_interface = data_interface;
            Xi = this.data_interface.Xi;
            n_r = size(Xi, 2);
            dist = zeros(n_r, n_r);
            for s = 1:n_r
                for k = 1:n_r
                    dist(s, k) = norm(Xi(:, s) - Xi(:, k))^2;
                end
            end
            this.C = exp(-dist);
            this.R = chol(this.C);
            this.Rinv = linsolve(this.R, eye(n_r));
        end

    end

end
