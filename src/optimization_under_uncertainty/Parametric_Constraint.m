%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Parametric_Constraint < Constraint
    % Define a parametric constraint function :math:`\c(\u, \z, \theta)\to\R^{n_u}` where
    % :math:`\u \in \R^{n_u}` is the state and
    % :math:`\z \in \R^{n_z}` is the control and
    % :math:`\theta \in \R^{n_theta}` are the parameters.

    properties
        theta_current                 % Current parameter vector
    end

    %% Constructor

    methods (Access = public)

        function this = Parametric_Constraint(theta)
            % Parameters
            % ----------
            % theta
            %   parameter, a vector.
            this@Constraint();
            this.theta_current = theta;
        end

    end

    %% Required abstract methods.

    methods (Abstract, Access = public)

        [u] = Parametric_State_Solve(this, z, theta)
        % Given :math:`\z` and `\theta`, solve the constraint equation
        % :math:`\c(\u,\z, \theta)=\0` for :math:`\u`, i.e., compute
        % :math:`\u = \S(\z, \theta)`.
        %
        % Parameters
        % ----------
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % theta
        %   Parameters :math:`\theta\in\R^{n_\theta}`.
        %
        % Returns
        % -------
        % u : vector
        %   State :math:`\u = \S(\z) \in \R^{n_u}`.

        [u_out] = Parametric_c_u_Transpose_Inverse_Apply(this, u_in, u, z, theta)
        % Compute the Jacobian-vector product :math:`c_u(\u,\z, \theta)\invtrp\v`,
        % i.e., solve :math:`c_u(\u,\z,\theta)\trp\bflambda = \v`
        % for :math:`\bflambda`.
        %
        % Parameters
        % ----------
        % u_in
        %   Search direction :math:`\v\in\R^{n_u}`.
        % u
        %   State :math:`\u\in\R^{n_u}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % theta
        %   Parameters :math:`\theta\in\R^{n_\theta}`.
        %
        % Returns
        % -------
        % u_out : vector
        %   Jacobian-vector product :math:`c_u(\u,\z,\theta)\invtrp\v\in\R^{n_u}`.

        [z_out] = Parametric_c_z_Transpose_Apply(this, u_in, u, z, theta)
        % Compute the Jacobian-vector product :math:`c_z(\u,\z,\theta)\trp\v`.
        %
        % Parameters
        % ----------
        % u_in
        %   Search direction :math:`\v\in\R^{n_u}`.
        % u
        %   State :math:`\u\in\R^{n_u}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % theta
        %   Parameters :math:`\theta\in\R^{n_\theta}`.
        %
        % Returns
        % -------
        % z_out : vector
        %   Jacobian-vector product :math:`c_z(\u,\z,\theta)\trp\v\in\R^{n_z}`.

        [u_out] = Parametric_c_u_Inverse_Apply(this, u_in, u, z, theta)
        % Compute the Jacobian-vector product :math:`c_u(\u,\z,\theta)^{-1}\v`,
        % i.e., solve :math:`c_u(\u,\z,\theta)\bfmu = \v` for :math:`\bfmu`.
        %
        % Parameters
        % ----------
        % u_in
        %   Search direction :math:`\v\in\R^{n_u}`.
        % u
        %   State :math:`\u\in\R^{n_u}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % theta
        %   Parameters :math:`\theta\in\R^{n_\theta}`.
        %
        % Returns
        % -------
        % u_out : vector
        %   Jacobian-vector product :math:`c_u(\u,\z,\theta)^{-1}\v\in\R^{n_u}`.

        [u_out] = Parametric_c_z_Apply(this, z_in, u, z, theta)
        % Compute the Jacobian-vector product :math:`c_z(\u,\z,\theta)\v`.
        %
        % Parameters
        % ----------
        % z_in
        %   Search direction :math:`\v\in\R^{n_z}`.
        % u
        %   State :math:`\u\in\R^{n_u}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % theta
        %   Parameters :math:`\theta\in\R^{n_\theta}`.
        %
        % Returns
        % -------
        % u_out : vector
        %   Jacobian-vector product :math:`c_z(\u,\z,\theta)\v\in\R^{n_u}`.

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
            % *Semi-abstract method.*
            % Explicitly form the constraint :math:`\c(\u,\z,\theta)`.
            % This method is only used for finite difference checks.
            %
            % Parameters
            % ----------
            % u
            %   State :math:`\u\in\R^{n_u}`.
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            % theta
            %   Parameters :math:`\theta\in\R^{n_\theta}`.
            %
            % Returns
            % -------
            % c : vector
            %   Constraint :math:`\c(\u,\z,\theta)\in\R^{n_u}`.
            con = error('Parametric_c() not implemented');
        end

        function [con] = c(this, u, z)
            con = this.Parametric_c(u, z, this.theta_current);
        end

        function [u_out] = Parametric_c_uu_Apply(this, u_in, u, z, lambda, theta)
            % *Semi-abstract method.*
            % Compute the vector-Hessian-vector product
            % :math:`\bflambda\trp c_{u,u}(\u,\z,\theta)\v`.
            %
            % Parameters
            % ----------
            % u_in
            %   Search direction :math:`\v\in\R^{n_u}`.
            % u
            %   State :math:`\u\in\R^{n_u}`.
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            % lambda
            %   Adjoint :math:`\bflambda\in\R^{n_u}`.
            % theta
            %   Parameters :math:`\theta\in\R^{n_\theta}`.
            %
            % Returns
            % -------
            % u_out : vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp c_{u,u}(\u,\z,\theta)\v\in\R^{n_u}`.
            u_out = error('Parametric_c_uu_Apply() not implemented');
        end

        function [u_out] = c_uu_Apply(this, u_in, u, z, lambda)
            u_out = this.Parametric_c_uu_Apply(u_in, u, z, lambda, this.theta_current);
        end

        function [u_out] = Parametric_c_uz_Apply(this, z_in, u, z, lambda, theta)
            % *Semi-abstract method.*
            % Compute the vector-Hessian-vector product
            % :math:`\bflambda\trp c_{u,z}(\u,\z,\theta)\v`.
            %
            % Parameters
            % ----------
            % z_in
            %   Search direction :math:`\v\in\R^{n_z}`.
            % u
            %   State :math:`\u\in\R^{n_u}`.
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            % lambda
            %   Adjoint :math:`\bflambda\in\R^{n_u}`.
            % theta
            %   Parameters :math:`\theta\in\R^{n_\theta}`.
            %
            % Returns
            % -------
            % u_out : vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp c_{u,z}(\u,\z,\theta)\v\in\R^{n_u}`.
            u_out = error('Parametric_c_uz_Apply() not implemented');
        end

        function [u_out] = c_uz_Apply(this, z_in, u, z, lambda)
            u_out = this.Parametric_c_uz_Apply(z_in, u, z, lambda, this.theta_current);
        end

        function [u_out] = Parametric_c_utheta_Apply(this, theta_in, u, z, lambda, theta)
            % *Semi-abstract method.*
            % Compute the vector-Hessian-vector product
            % :math:`\bflambda\trp c_{u,\theta}(\u,\z,\theta)\v`.
            %
            % Parameters
            % ----------
            % theta_in
            %   Search direction :math:`\v\in\R^{n_\theta}`.
            % u
            %   State :math:`\u\in\R^{n_u}`.
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            % lambda
            %   Adjoint :math:`\bflambda\in\R^{n_u}`.
            % theta
            %   Parameters :math:`\theta\in\R^{n_\theta}`.
            %
            % Returns
            % -------
            % u_out : vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp c_{u,\theta}(\u,\z,\theta)\v\in\R^{n_u}`.
            u_out = error('Parametric_c_utheta_Apply() not implemented');
        end

        function [z_out] = Parametric_c_zu_Apply(this, u_in, u, z, lambda, theta)
            % *Semi-abstract method.*
            % Compute the vector-Hessian-vector product
            % :math:`\bflambda\trp c_{z,u}(\u,\z,\theta)\v`.
            %
            % Parameters
            % ----------
            % u_in
            %   Search direction :math:`\v\in\R^{n_u}`.
            % u
            %   State :math:`\u\in\R^{n_u}`.
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            % lambda
            %   Adjoint :math:`\bflambda\in\R^{n_u}`.
            % theta
            %   Parameters :math:`\theta\in\R^{n_\theta}`.
            %
            % Returns
            % -------
            % z_out : vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp c_{z,u}(\u,\z,\theta)\v\in\R^{n_z}`.
            z_out = error('Parametric_c_zu_Apply() not implemented');
        end

        function [z_out] = c_zu_Apply(this, u_in, u, z, lambda)
            z_out = this.Parametric_c_zu_Apply(u_in, u, z, lambda, this.theta_current);
        end

        function [z_out] = Parametric_c_zz_Apply(this, z_in, u, z, lambda, theta)
            % *Semi-abstract method.*
            % Compute the vector-Hessian-vector product
            % :math:`\bflambda\trp c_{z,z}(\u,\z,\theta)\v`.
            %
            % Parameters
            % ----------
            % z_in
            %   Search direction :math:`\v\in\R^{n_z}`.
            % u
            %   State :math:`\u\in\R^{n_u}`.
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            % lambda
            %   Adjoint :math:`\bflambda\in\R^{n_u}`.
            % theta
            %   Parameters :math:`\theta\in\R^{n_\theta}`.
            %
            % Returns
            % -------
            % z_out : vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp c_{z,z}(\u,\z,\theta)\v\in\R^{n_z}`.
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
            % Check the implementation of the following via finite differences.
            %
            % * :meth:`c_z_Apply()` for :math:`\c_z(\u,\z)\v`.
            % * :meth:`c_u_Inverse_Apply()` for :math:`\c_u(\u,\z)^{-1}\v`.
            %
            % Also check that the following functions are consistent.
            %
            % * :meth:`c_z_Apply()` and :meth:`c_z_Transpose_Apply()`.
            % * :meth:`c_u_Inverse_Apply()` and :meth:`c_u_Transpose_Inverse_Apply()`.
            %
            % Finally, check that :meth:`State_Solve()` and :meth:`c()` are inverses.
            %
            % Note
            % ----
            % This check requires :meth:`c()` to be implemented.
            %
            % Parameters
            % ----------
            % u
            %   State :math:`\mathbf{u} \in \mathbb{R}^{n_u}`.
            % z
            %   Control :math:`\mathbf{z} \in \mathbb{R}^{n_z}`.
            %
            % Returns
            % -------
            % diffs_z : vector
            %   Finite difference errors for :math:`\c_z(\u, \z)\v`.
            % jacobian_z_transpose_check : double
            %   Error from comparing :math:`\c_z(\u,\z)\trp\v` with :math:`\c_z(\u,\z)\v`.
            % diffs_u : vector
            %   Finite difference errors for :math:`\c_u(\u, \z)\v`.
            % jacobian_u_transpose_check : double
            %   Error from comparing :math:`\c_u(\u,\z)\invtrp\v` with :math:`\c_u(\u,\z)^{-1}\v`.
            % solve_res : double
            %   Error from comparing :math:`\S(\z)` with :math:`\c(\u, \z)`.

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
