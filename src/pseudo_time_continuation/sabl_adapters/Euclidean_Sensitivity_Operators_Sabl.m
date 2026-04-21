%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Euclidean_Sensitivity_Operators_Sabl < Euclidean_Sensitivity_Operators

    properties
        obj                 % Instance of a subclass of :class:`Objective`.
        pcon                % Instance of a subclass of :class:`Parameterized_Constraint`.
        current_u           % Current state.
        current_z           % Current control.
        current_theta       % Current parameters.
        current_lambda      % Current adjoint.
        verbose             % Verbosity.
        Gauss_Newton_Hess   % Use Gauss-Newton approximation of the Hessian.
    end

    methods

        function this = Euclidean_Sensitivity_Operators_Sabl(obj, pcon)
            % Parameters
            % ----------
            % obj
            %   Objective function, an instance of a subclass of :class:`Objective`.
            % pcon
            %   Constraint equations, an instance of a subclass of :class:`Parameterized_Constraint`.
            this.obj = obj;
            this.pcon = pcon;
            this.verbose = true;
            this.Gauss_Newton_Hess = false;
        end

        function [] = Update(this, z, theta)
            if norm(z - this.current_z) ~= 0 || norm(theta - this.current_theta) ~= 0
                this.Gradient(z, theta);
            end
        end

        function [grad, val] = Euclidean_Gradient(this, z, theta)
            if isempty(this.current_z)
                u = this.pcon.Parameterized_State_Solve(z, theta);
                [val, grad_u, grad_z] = this.obj.J(u, z);
                lambda = this.pcon.Parameterized_c_u_Transpose_Inverse_Apply(-grad_u, u, z, theta);
                grad = this.pcon.Parameterized_c_z_Transpose_Apply(lambda, u, z, theta);
                grad = grad + grad_z;
                this.current_u = u;
                this.current_z = z;
                this.current_theta = theta;
                this.current_lambda = lambda;
            elseif norm(z - this.current_z) ~= 0 || norm(theta - this.current_theta) ~= 0
                u = this.pcon.Parameterized_State_Solve(z, theta);
                [val, grad_u, grad_z] = this.obj.J(u, z);
                lambda = this.pcon.Parameterized_c_u_Transpose_Inverse_Apply(-grad_u, u, z, theta);
                grad = this.pcon.Parameterized_c_z_Transpose_Apply(lambda, u, z, theta);
                grad = grad + grad_z;
                this.current_u = u;
                this.current_z = z;
                this.current_theta = theta;
                this.current_lambda = lambda;
            else
                u = this.current_u;
                [val, ~, grad_z] = this.obj.J(u, z);
                lambda = this.current_lambda;
                grad = this.pcon.Parameterized_c_z_Transpose_Apply(lambda, u, z, theta);
                grad = grad + grad_z;
            end
        end

        function [z_out] = Euclidean_Apply_Hessian(this, z_in, z, theta)
            this.Update(z, theta);
            w = this.pcon.Parameterized_c_z_Apply(z_in, this.current_u, this.current_z, this.current_theta);
            mu = this.pcon.Parameterized_c_u_Inverse_Apply(-w, this.current_u, this.current_z, this.current_theta);
            yJ = this.obj.J_uu_Apply(mu, this.current_u, this.current_z) + this.obj.J_uz_Apply(z_in, this.current_u, this.current_z);
            xJ = this.obj.J_zu_Apply(mu, this.current_u, this.current_z) + this.obj.J_zz_Apply(z_in, this.current_u, this.current_z);
            if this.Gauss_Newton_Hess
                gamma = this.pcon.Parameterized_c_u_Transpose_Inverse_Apply(-yJ, this.current_u, this.current_z, this.current_theta);
                xc = this.pcon.Parameterized_c_z_Transpose_Apply(gamma, this.current_u, this.current_z, this.current_theta);
            else
                yc = this.pcon.Parameterized_c_uu_Apply(mu, this.current_u, this.current_z, this.current_lambda, this.current_theta) + this.pcon.Parameterized_c_uz_Apply(z_in, this.current_u, this.current_z, this.current_lambda, this.current_theta);
                gamma = this.pcon.Parameterized_c_u_Transpose_Inverse_Apply(-(yJ + yc), this.current_u, this.current_z, this.current_theta);
                xc = this.pcon.Parameterized_c_z_Transpose_Apply(gamma, this.current_u, this.current_z, this.current_theta);
                xc = xc + this.pcon.Parameterized_c_zu_Apply(mu, this.current_u, this.current_z, this.current_lambda, this.current_theta);
                xc = xc + this.pcon.Parameterized_c_zz_Apply(z_in, this.current_u, this.current_z, this.current_lambda, this.current_theta);
            end
            z_out = xJ + xc;
        end

        function [z_out] = Euclidean_Apply_B(this, theta_in, z, theta)
            this.Update(z, theta);
            w = this.pcon.Parameterized_c_theta_Apply(theta_in, this.current_u, this.current_z, this.current_theta);
            xi = this.pcon.Parameterized_c_u_Inverse_Apply(-w, this.current_u, this.current_z, this.current_theta);
            yJ = this.obj.J_uu_Apply(xi, this.current_u, this.current_z);
            xJ = this.obj.J_zu_Apply(xi, this.current_u, this.current_z);
            yc = this.pcon.Parameterized_c_uu_Apply(xi, this.current_u, this.current_z, this.current_lambda, this.current_theta) + this.pcon.Parameterized_c_utheta_Apply(theta_in, this.current_u, this.current_z, this.current_lambda, this.current_theta);
            beta = this.pcon.Parameterized_c_u_Transpose_Inverse_Apply(-(yJ + yc), this.current_u, this.current_z, this.current_theta);
            xc = this.pcon.Parameterized_c_z_Transpose_Apply(beta, this.current_u, this.current_z, this.current_theta);
            xc = xc + this.pcon.Parameterized_c_zu_Apply(xi, this.current_u, this.current_z, this.current_lambda, this.current_theta);
            xc = xc + this.pcon.Parameterized_c_ztheta_Apply(theta_in, this.current_u, this.current_z, this.current_lambda, this.current_theta);
            z_out = xJ + xc;
        end

    end
end
