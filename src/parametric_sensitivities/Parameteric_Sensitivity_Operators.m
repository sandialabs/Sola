classdef Parameteric_Sensitivity_Operators < handle

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

        function this = Parameteric_Sensitivity_Operators(obj, pcon)
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

        function [grad] = Solve_Forward_and_Adjoint_Problems(this, z, theta)
            u = this.pcon.Parameterized_State_Solve(z, theta);
            [~, grad_u, grad_z] = this.obj.J(u, z);
            lambda = this.pcon.Parameterized_c_u_Transpose_Inverse_Apply(-grad_u, u, z, theta);
            grad = this.pcon.Parameterized_c_z_Transpose_Apply(lambda, u, z, theta);
            grad = grad + grad_z;
            this.current_u = u;
            this.current_z = z;
            this.current_theta = theta;
            this.current_lambda = lambda;
        end

        function [z_out] = Apply_RS_Hessian(this, z_in)
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

        function [z_out] = Apply_B(this, theta_in)
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

        function [diffs] = Finite_Difference_Hessian_Check(this, z, theta)
            grad = this.Solve_Forward_and_Adjoint_Problems(z, theta);
            n = length(grad);
            v = randn(n, 1);
            v = v / norm(v);
            Hv = this.Apply_RS_Hessian(v);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_hv = zeros(n, p);
            diffs = zeros(p, 1);
            for k = 1:p
                gradk = this.Solve_Forward_and_Adjoint_Problems(z + h(k) * v, theta);
                fd_hv(:, k) = (gradk - grad) / h(k);
                diffs(k) = norm(fd_hv(:, k) - Hv) / norm(Hv);
            end
            if this.verbose
                disp('Hessian finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

        function [diffs] = Finite_Difference_B_Check(this, z, theta)
            grad = this.Solve_Forward_and_Adjoint_Problems(z, theta);
            n = length(theta);
            v = randn(n, 1);
            v = v / norm(v);
            Bv = this.Apply_B(v);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_Bv = zeros(length(z), p);
            diffs = zeros(p, 1);
            for k = 1:p
                gradk = this.Solve_Forward_and_Adjoint_Problems(z, theta + h(k) * v);
                fd_Bv(:, k) = (gradk - grad) / h(k);
                diffs(k) = norm(fd_Bv(:, k) - Bv) / norm(Bv);
            end
            if this.verbose
                disp('B finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

    end
end
