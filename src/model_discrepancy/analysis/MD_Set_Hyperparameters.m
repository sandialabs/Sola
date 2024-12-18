classdef MD_Set_Hyperparameters < handle

    properties
        u_prior_interface
        z_prior_interface
        data_interface
    end

    methods

        function [alpha_d] = Determine_alpha_d(this)
            alpha_d = (0.001 * mean(abs(this.data_interface.D(:))))^2;
        end

        function [alpha_z] = Determine_alpha_z(this, z)
            E_z = this.z_prior_interface.Apply_E_z_Inverse(z - this.data_interface.z_opt);
            M_E_z = this.z_prior_interface.Apply_M_z(E_z);
            alpha_z = (0.25^2) / (E_z' * M_E_z);
        end

        function [alpha_u] = Determine_alpha_u(this)
            delta = this.data_interface.Load_d_Data();
            delta_norm = sqrt(delta' * this.u_prior_interface.Apply_M_u(delta));
            delta_norm = mean(diag(delta_norm));

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

        function this = MD_Set_Hyperparameters(u_prior_interface, z_prior_interface, data_interface)
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.data_interface = data_interface;
        end

    end

end
