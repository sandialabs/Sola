%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Parametric_Constraint < Constraint
    % Define a parametric constraint function c(u,z,theta)

    properties
        theta_current                 % Current parameter vector
    end

    %% Constructor

    methods (Access = public)

        function this = Parametric_Constraint(theta)
            arguments
                theta (:, 1) {mustBeNumeric}
            end
            this@Constraint();
            this.theta_current = theta;
        end

    end

    %% Required abstract methods.

    methods (Abstract, Access = public)

        [u] = Parametric_State_Solve(this, z, theta)

        [u_out] = Parametric_c_u_Transpose_Inverse_Apply(this, u_in, u, z, theta)

        [z_out] = Parametric_c_z_Transpose_Apply(this, u_in, u, z, theta)

        [u_out] = Parametric_c_u_Inverse_Apply(this, u_in, u, z, theta)

        [u_out] = Parametric_c_z_Apply(this, z_in, u, z, theta)

    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u = this.Parametric_State_Solve(z, this.theta_current);
        end

        function [u_out] = c_u_Transpose_Inverse_Apply(this, u_in, u, z)
            u_out = this.Parametric_c_u_Transpose_Inverse_Apply(u_in, u, z, this.theta_current);
        end

        function [z_out] = c_z_Transpose_Apply(this, u_in, u, z)
            z_out = this.Parametric_c_z_Transpose_Apply(u_in, u, z, this.theta_current);
        end

        function [u_out] = c_u_Inverse_Apply(this, u_in, u, z)
            u_out = this.Parametric_c_u_Inverse_Apply(u_in, u, z, this.theta_current);
        end

        function [u_out] = c_z_Apply(this, z_in, u, z)
            u_out = this.Parametric_c_z_Apply(z_in, u, z, this.theta_current);
        end

    end

    %% Semi-abstract methods, required when Gauss_Newton_Hess = false.

    methods (Access = public)

        function [con] = Parametric_c(this, u, z, theta)
            con = error('Parametric_c() not implemented');
        end

        function [con] = c(this, u, z)
            con = this.Parametric_c(u, z, this.theta_current);
        end

        function [u_out] = Parametric_c_uu_Apply(this, u_in, u, z, lambda, theta)
            u_out = error('Parametric_c_uu_Apply() not implemented');
        end

        function [u_out] = c_uu_Apply(this, u_in, u, z, lambda)
            u_out = this.Parametric_c_uu_Apply(u_in, u, z, lambda, this.theta_current);
        end

        function [u_out] = Parametric_c_uz_Apply(this, z_in, u, z, lambda, theta)
            u_out = error('Parametric_c_uz_Apply() not implemented');
        end

        function [u_out] = c_uz_Apply(this, z_in, u, z, lambda)
            u_out = this.Parametric_c_uz_Apply(z_in, u, z, lambda, this.theta_current);
        end

        function [u_out] = Parametric_c_utheta_Apply(this, theta_in, u, z, lambda, theta)
            u_out = error('Parametric_c_utheta_Apply() not implemented');
        end

        function [z_out] = Parametric_c_zu_Apply(this, u_in, u, z, lambda, theta)
            z_out = error('Parametric_c_zu_Apply() not implemented');
        end

        function [z_out] = c_zu_Apply(this, u_in, u, z, lambda)
            z_out = this.Parametric_c_zu_Apply(u_in, u, z, lambda, this.theta_current);
        end

        function [z_out] = Parametric_c_zz_Apply(this, z_in, u, z, lambda, theta)
            z_out = error('Parametric_c_zz_Apply() not implemented');
        end

        function [z_out] = c_zz_Apply(this, z_in, u, z, lambda)
            z_out = this.Parametric_c_zz_Apply(z_in, u, z, lambda, this.theta_current);
        end

    end

    %% Finite difference checks.

    methods (Access = public)

        function [diffs_z, jacobian_z_transpose_check, diffs_u, jacobian_u_transpose_check, solve_res] = Parametric_Finite_Difference_Constraint_Check(this, u, z, theta)
            theta_tmp = this.theta_current;
            this.theta_current = theta;
            [diffs_z, jacobian_z_transpose_check, diffs_u, jacobian_u_transpose_check, solve_res] = this.Finite_Difference_Constraint_Check(u, z);
            this.theta_current = theta_tmp;
        end

        function [diffs_z, jacobian_z_transpose_check, diffs_u, jacobian_u_transpose_check, solve_res] = Finite_Difference_Constraint_Check(this, u, z)

            % c_z_Apply check
            c = this.c(u, z);
            h = 10.^(-2:-1:-6);
            p = length(h);
            v = randn(length(z), 1);
            cz_v = this.c_z_Apply(v, u, z);
            diffs_z = zeros(p, 1);
            for k = 1:p
                ck =  this.c(u, z + h(k) * v);
                diffs_z(k) = norm(cz_v - (ck - c) / h(k)) / norm(cz_v);
            end
            disp('Constraint z Jacobian finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k), '%.1e'), ' and error = ', num2str(diffs_z(k))]);
            end
            disp(' ');

            % c_z_Transpose_Apply check
            vu = randn(length(u), 1);
            vu_cz_v = vu' * cz_v;
            ctz_vu = this.c_z_Transpose_Apply(vu, u, z);
            v_ctz_vu = v' * ctz_vu;
            jacobian_z_transpose_check = abs(vu_cz_v - v_ctz_vu) / abs(vu_cz_v);
            if jacobian_z_transpose_check > 1.e-13
                disp('Error in c_z_Transpose_Apply');
            end

            % c_u_Inverse_Apply check
            diffs_u = zeros(p, 1);
            for k = 1:p
                cu_vu_fd = (this.c(u + h(k) * vu, z) - c) / h(k);
                diffs_u(k) = norm(this.c_u_Inverse_Apply(cu_vu_fd, u, z) - vu) / norm(vu);
            end
            disp('Constraint u Jacobian Inverse finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_u(k))]);
            end
            disp(' ');

            % c_u_Transpose_Inverse_Apply check
            cu_inv_vu = this.c_u_Inverse_Apply(vu, u, z);
            d = randn(length(u), 1);
            cu_trans_inv_d = this.c_u_Transpose_Inverse_Apply(d, u, z);
            jacobian_u_transpose_check = abs(vu' * cu_trans_inv_d - d' * cu_inv_vu) / abs(d' * cu_inv_vu);
            if jacobian_u_transpose_check > 1.e-13
                disp('Error in c_u_Transpose_Inverse_Apply');
            end

            % State solve check
            u_solve = this.State_Solve(z);
            solve_res = norm(this.c(u_solve, z));
            if solve_res > 1.e-10
                disp('Error in state solve');
            end
        end

    end
end
