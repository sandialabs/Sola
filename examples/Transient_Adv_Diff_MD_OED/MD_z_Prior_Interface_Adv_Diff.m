classdef MD_z_Prior_Interface_Adv_Diff < MD_z_Prior_Interface

    properties
        Wz
        F
        alpha_z
    end

    methods (Access = public)

        function [z_out] = Apply_W_z_Inverse(this, z_in)
            z_out = this.alpha_z * this.F * (this.F' * z_in);
        end

        % Compute samples from a mean zero Gaussian with covariance W_z^{-1}
        function [z_out] = Sample_with_Covariance_W_z_Inverse(this, num_samples)
            z_out = sqrt(this.alpha_z) * this.F * randn(size(this.F, 1), num_samples);
        end

        % Apply W_z matrix
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_W_z(this, z_in)
            z_out = (1 / this.alpha_z) * this.Wz * z_in;
        end

        function this = MD_z_Prior_Interface_Adv_Diff(obj)
            this.Wz = kron(diag(obj.time_weights(2:end)), obj.Br' * obj.M * obj.Br);
            R = chol(this.Wz);
            this.F = linsolve(R, eye(size(R, 1)));

            this.alpha_z = 5.e-5;
        end

    end

end
