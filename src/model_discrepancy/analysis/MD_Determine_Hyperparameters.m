classdef MD_Determine_Hyperparameters < handle

    properties
        u_prior_interface
        z_prior_interface
        data_interface
    end

    methods

        function [alpha_d] = Determine_alpha_d(this, scaling)
            arguments
                this
                scaling(1,1) double {mustBeFinite} = 0.001;
            end
            alpha_d = (scaling * mean(abs(this.data_interface.D(:))))^2;
        end

        function [alpha_z] = Determine_alpha_z(this, scaling)
            arguments
                this
                scaling(1,1) double {mustBeFinite} = 0.2;
            end

            num_samples = 50;
            samples = this.z_prior_interface.Sample_with_Covariance_W_z_Inverse(num_samples);
            tmp = this.z_prior_interface.Apply_M_z(this.data_interface.z_opt);
            tmp = sqrt(tmp'*this.data_interface.z_opt);
            for k = 1:num_samples
                samples(:,k) = tmp*samples(:,k)/sqrt(samples(:,k)'*this.z_prior_interface.Apply_M_z(samples(:,k)));
            end
            sample_norms = zeros(num_samples,1);
            for k = 1:num_samples
                E_z = this.z_prior_interface.Apply_E_z_Inverse(samples(:,k));
                sample_norms(k) = E_z'*this.z_prior_interface.Apply_M_z(E_z);
            end
            alpha_z = (scaling^2) / mean(sample_norms);
        end

        function [alpha_u] = Determine_alpha_u(this, scaling)
            arguments
                this
                scaling(1,1) double {mustBeFinite} = 1.0;
            end
            delta = this.data_interface.Load_d_Data();
            delta_norm = sqrt(delta' * this.u_prior_interface.Apply_M_u(delta));
            delta_norm = scaling * mean(diag(delta_norm));

            target = delta_norm;
            dof = sum(this.u_prior_interface.sing_vals.^2);

            % Want to solve:
            % sqrt(2) * Gamma(0.5 * (alpha_u * dof + 1) ) = target * Gamma(0.5 * alpha_u * dof)
            fun = @(alpha_u) this.alpha_u_fun(alpha_u, dof, target);
            [alpha_u, ~, flag] = fsolve(fun, this.u_prior_interface.alpha_u, optimoptions('fsolve', 'Display', 'off'));
            if flag <= 0
                disp('Error in nonlinear solve to determine alpha_u');
            end
        end

        function [val] = alpha_u_fun(this, alpha_u, dof, target)
            val = sqrt(2) * gamma(0.5 * (alpha_u * dof + 1)) - target * gamma(0.5 * alpha_u * dof);
        end

        function this = MD_Determine_Hyperparameters(u_prior_interface, z_prior_interface, data_interface)
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.data_interface = data_interface;
        end

    end

end
