%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Parameterized_Constraint < Constraint
    % Define a parameterized constraint function :math:`\c(\u, \z, \theta)\to\R^{n_u}` where
    % :math:`\u \in \R^{n_u}` is the state and
    % :math:`\z \in \R^{n_z}` is the control and
    % :math:`\theta \in \R^{n_theta}` are the parameters.

    properties
        theta_current                 % Current parameter vector
    end

    %% Constructor

    methods (Access = public)

        function this = Parameterized_Constraint(theta)
            % Parameters
            % ----------
            % theta
            %   parameter, a vector.
            this.theta_current = theta;
        end

    end

    %% Required abstract methods.

    methods (Abstract, Access = public)

        [u] = Parameterized_State_Solve(this, z, theta)
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

        [u_out] = Parameterized_c_u_Transpose_Inverse_Apply(this, u_in, u, z, theta)
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

        [z_out] = Parameterized_c_z_Transpose_Apply(this, u_in, u, z, theta)
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

        [u_out] = Parameterized_c_u_Inverse_Apply(this, u_in, u, z, theta)
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

        [u_out] = Parameterized_c_z_Apply(this, z_in, u, z, theta)
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

        [u_out] = Parameterized_c_theta_Apply(this, theta_in, u, z, theta)
        % Compute the Jacobian-vector product :math:`c_z(\u,\z,\theta)\v`.
        %
        % Parameters
        % ----------
        % theta_in
        %   Search direction :math:`\v\in\R^{n_\theta}`.
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
        %   Jacobian-vector product :math:`c_\theta(\u,\z,\theta)\v\in\R^{n_u}`.

    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u = this.Parameterized_State_Solve(z, this.theta_current);
        end

        function [u_out] = c_u_Transpose_Inverse_Apply(this, u_in, u, z)
            u_out = this.Parameterized_c_u_Transpose_Inverse_Apply(u_in, u, z, this.theta_current);
        end

        function [z_out] = c_z_Transpose_Apply(this, u_in, u, z)
            z_out = this.Parameterized_c_z_Transpose_Apply(u_in, u, z, this.theta_current);
        end

        function [u_out] = c_u_Inverse_Apply(this, u_in, u, z)
            u_out = this.Parameterized_c_u_Inverse_Apply(u_in, u, z, this.theta_current);
        end

        function [u_out] = c_z_Apply(this, z_in, u, z)
            u_out = this.Parameterized_c_z_Apply(z_in, u, z, this.theta_current);
        end

    end

    %% Semi-abstract methods, required when Gauss_Newton_Hess = false.

    methods (Access = public)

        function [con] = Parameterized_c(this, u, z, theta)
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
            con = error('Parameterized_c() not implemented');
        end

        function [con] = c(this, u, z)
            con = this.Parameterized_c(u, z, this.theta_current);
        end

        function [u_out] = Parameterized_c_uu_Apply(this, u_in, u, z, lambda, theta)
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
            u_out = error('Parameterized_c_uu_Apply() not implemented');
        end

        function [u_out] = c_uu_Apply(this, u_in, u, z, lambda)
            u_out = this.Parameterized_c_uu_Apply(u_in, u, z, lambda, this.theta_current);
        end

        function [u_out] = Parameterized_c_uz_Apply(this, z_in, u, z, lambda, theta)
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
            u_out = error('Parameterized_c_uz_Apply() not implemented');
        end

        function [u_out] = c_uz_Apply(this, z_in, u, z, lambda)
            u_out = this.Parameterized_c_uz_Apply(z_in, u, z, lambda, this.theta_current);
        end

        function [u_out] = Parameterized_c_utheta_Apply(this, theta_in, u, z, lambda, theta)
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
            u_out = error('Parameterized_c_utheta_Apply() not implemented');
        end

        function [z_out] = Parameterized_c_zu_Apply(this, u_in, u, z, lambda, theta)
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
            z_out = error('Parameterized_c_zu_Apply() not implemented');
        end

        function [z_out] = c_zu_Apply(this, u_in, u, z, lambda)
            z_out = this.Parameterized_c_zu_Apply(u_in, u, z, lambda, this.theta_current);
        end

        function [z_out] = Parameterized_c_ztheta_Apply(this, theta_in, u, z, lambda, theta)
            % *Semi-abstract method.*
            % Compute the vector-Hessian-vector product
            % :math:`\bflambda\trp c_{z,\theta}(\u,\z,\theta)\v`.
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
            % z_out : vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp c_{z,\theta}(\u,\z,\theta)\v\in\R^{n_z}`.
            z_out = error('Parameterized_c_ztheta_Apply() not implemented');
        end

        function [z_out] = Parameterized_c_zz_Apply(this, z_in, u, z, lambda, theta)
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
            z_out = error('Parameterized_c_zz_Apply() not implemented');
        end

        function [z_out] = c_zz_Apply(this, z_in, u, z, lambda)
            z_out = this.Parameterized_c_zz_Apply(z_in, u, z, lambda, this.theta_current);
        end

    end

    %% Finite difference checks.

    methods (Access = public)

        function [diffs_z, jacobian_z_transpose_check, diffs_u, jacobian_u_transpose_check, diffs_theta, solve_res] = Parameterized_Finite_Difference_Constraint_Check(this, u, z, theta)
            theta_tmp = this.theta_current;
            this.theta_current = theta;
            [diffs_z, jacobian_z_transpose_check, diffs_u, jacobian_u_transpose_check, solve_res] = this.Finite_Difference_Constraint_Check(u, z);
            this.theta_current = theta_tmp;

            % Need to implement theta derivative checks

            % Parameterized_c_theta_Apply check
            c = this.Parameterized_c(u, z, theta);
            h = 10.^(-2:-1:-6);
            p = length(h);
            v = randn(length(theta), 1);
            ctheta_v = this.Parameterized_c_theta_Apply(v, u, z);
            diffs_theta = zeros(p, 1);
            for k = 1:p
                ck =  this.Parameterized_c(u, z, theta + h(k) * v);
                diffs_theta(k) = norm(ctheta_v - (ck - c) / h(k)) / norm(ctheta_v);
            end
            disp('Constraint theta Jacobian finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k), '%.1e'), ' and error = ', num2str(diffs_theta(k))]);
            end
            disp(' ');

        end

    end
end
