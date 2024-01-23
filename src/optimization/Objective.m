classdef Objective < handle
    % Define a scalar-valued objective function
    %
    % .. math:: J(\u, \z) \to \R
    %
    % where
    % :math:`\u \in \R^{n_u}` is the state and
    % :math:`\z \in \R^{n_z}` is the control.

    methods (Abstract, Access = public)

        [val, grad_u, grad_z] = J(this, u, z)
        % Evaluate the objective function :math:`J(\u,\z)` and its
        % gradients :math:`\grad{u}J(\u,\z)` and :math:`\grad{z}J(\u,\z)`.
        %
        % Parameters
        % ----------
        % u
        %   State :math:`\u\in\R^{n_u}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        %
        % Returns
        % -------
        % val : double
        %   Objective value :math:`J(\u,\z)`.
        % grad_u : vector
        %   Objective gradient :math:`\grad{u}J(\u,\z)`.
        % grad_z : vector
        %   Objective gradient :math:`\grad{z}J(\u,\z)`.

        [u_out] = J_uu_Apply(this, u_in, u, z)
        % Compute :math:`\grad{u,u}J(\u,\z)\v`.
        %
        % Parameters
        % ----------
        % u_in
        %   Search direction :math:`\v\in\R^{n_u}`.
        % u
        %   State :math:`\u\in\R^{n_u}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        %
        % Returns
        % -------
        % u_out : vector
        %   Gradient-vector product :math:`\grad{u,u}J(\u,\z)\v\in\R^{n_u}`.

        [u_out] = J_uz_Apply(this, z_in, u, z)
        % Compute :math:`\grad{u,z}J(\u,\z)\v`.
        %
        % Parameters
        % ----------
        % z_in
        %   Search direction :math:`\v\in\R^{n_z}`.
        % u
        %   State :math:`\u\in\R^{n_u}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        %
        % Returns
        % -------
        % u_out : vector
        %   Gradient-vector product :math:`\grad{u,z}J(\u,\z)\v\in\R^{n_u}`.

        [z_out] = J_zu_Apply(this, u_in, u, z)
        % Compute :math:`\grad{z,u}J(\u,\z)\v`.
        %
        % Parameters
        % ----------
        % u_in
        %   Search direction :math:`\v\in\R^{n_u}`.
        % u
        %   State :math:`\u\in\R^{n_u}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        %
        % Returns
        % -------
        % z_out : vector
        %   Gradient-vector product :math:`\grad{z,u}J(\u,\z)\v\in\R^{n_z}`.

        [z_out] = J_zz_Apply(this, z_in, u, z)
        % Compute :math:`\grad{z,z}J(\u,\z)\v`.
        %
        % Parameters
        % ----------
        % z_in
        %   Search direction :math:`\v\in\R^{n_z}`.
        % u
        %   State :math:`\u\in\R^{n_u}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        %
        % Returns
        % -------
        % z_out : vector
        %   Gradient-vector product :math:`\grad{z,z}J(\u,\z)\v\in\R^{n_z}`.

    end

    methods (Access = public)

        function this = Objective()

        end

        % Finite difference test functions

        function [diffs_u, diffs_z] = Finite_Difference_Gradient_Check(this, u, z)
            % Check the implementation of :meth:`J()` via finite differences.
            %
            % Parameters
            % ----------
            % u
            %   State :math:`\u \in \R^{n_u}`.
            % z
            %   Control :math:`\z \in \R^{n_z}`.
            %
            % Returns
            % -------
            % diffs_u : vector
            %   Finite difference errors for :math:`\grad{u}J(\u,\z)`.
            % diffs_z : vector
            %   Finite difference errors for :math:`\grad{z}J(\u,\z)`.
            [val, grad_u, grad_z] = this.J(u, z);

            h = 10.^(-2:-1:-6);
            p = length(h);

            m = length(grad_u);
            du = randn(m, 1);
            du = du / norm(du);
            grad_du = du' * grad_u;
            fd_grad = zeros(p, 1);
            for k = 1:p
                valk = this.J(u + h(k) * du, z);
                fd_grad(k) = (valk - val) / h(k);
            end

            diffs_u = abs(grad_du - fd_grad) / abs(grad_du);
            disp('u gradient finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_u(k))]);
            end
            disp(' ');

            n = length(grad_z);
            dz = randn(n, 1);
            dz = dz / norm(dz);
            grad_dz = dz' * grad_z;
            fd_grad = zeros(p, 1);
            for k = 1:p
                valk = this.J(u, z + h(k) * dz);
                fd_grad(k) = (valk - val) / h(k);
            end

            diffs_z = abs(grad_dz - fd_grad) / abs(grad_dz);
            disp('z gradient finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_z(k))]);
            end
            disp(' ');

        end

        function [diffs_uu, diffs_uz, diffs_zu, diffs_zz] = Finite_Difference_Hessian_Check(this, u, z)
            % Check the implementation of the following via finite differences.
            %
            % * :meth:`J_uu_Apply()` for :math:`\grad{u,u}J(\u,\z)`.
            % * :meth:`J_uz_Apply()` for :math:`\grad{u,z}J(\u,\z)`.
            % * :meth:`J_zu_Apply()` for :math:`\grad{z,u}J(\u,\z)`.
            % * :meth:`J_zz_Apply()` for :math:`\grad{z,z}J(\u,\z)`.
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
            % diffs_uu : vector
            %   Finite difference errors for :math:`\grad{u,u}J(\u,\z)`.
            % diffs_uz : vector
            %   Finite difference errors for :math:`\grad{u,z}J(\u,\z)`.
            % diffs_zu : vector
            %   Finite difference errors for :math:`\grad{z,u}J(\u,\z)`.
            % diffs_zz : vector
            %   Finite difference errors for :math:`\grad{z,z}J(\u,\z)`.

            [~, grad_u, grad_z] = this.J(u, z);
            h = 10.^(-2:-1:-6);
            p = length(h);
            m = length(grad_u);
            n = length(grad_z);

            % Check J_uu_Apply().
            v = randn(m, 1);
            v = v / norm(v);
            Hv = this.J_uu_Apply(v, u, z);
            fd_hv = zeros(m, p);
            diffs_uu = zeros(p, 1);
            for k = 1:p
                [~, gradk] = this.J(u + h(k) * v, z);
                fd_hv(:, k) = (gradk - grad_u) / h(k);
                diffs_uu(k) = norm(fd_hv(:, k) - Hv);
            end
            if norm(Hv) > 0
                diffs_uu = diffs_uu / norm(Hv);
            end
            disp('uu Hessian finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_uu(k))]);
            end
            disp(' ');

            % Check J_uz_Apply().
            v = randn(n, 1);
            v = v / norm(v);
            Hv = this.J_uz_Apply(v, u, z);
            fd_hv = zeros(m, p);
            diffs_uz = zeros(p, 1);
            for k = 1:p
                [~, gradk] = this.J(u, z + h(k) * v);
                fd_hv(:, k) = (gradk - grad_u) / h(k);
                diffs_uz(k) = norm(fd_hv(:, k) - Hv);
            end
            if norm(Hv) > 0
                diffs_uz = diffs_uz / norm(Hv);
            end
            disp('uz Hessian finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_uz(k))]);
            end
            disp(' ');

            % Check J_zu_Apply().
            v = randn(m, 1);
            v = v / norm(v);
            Hv = this.J_zu_Apply(v, u, z);
            fd_hv = zeros(n, p);
            diffs_zu = zeros(p, 1);
            for k = 1:p
                [~, ~, gradk] = this.J(u + h(k) * v, z);
                fd_hv(:, k) = (gradk - grad_z) / h(k);
                diffs_zu(k) = norm(fd_hv(:, k) - Hv);
            end
            if norm(Hv) > 0
                diffs_zu = diffs_zu / norm(Hv);
            end
            disp('zu Hessian finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_zu(k))]);
            end
            disp(' ');

            % Check J_zz_Apply().
            v = randn(n, 1);
            v = v / norm(v);
            Hv = this.J_zz_Apply(v, u, z);
            fd_hv = zeros(n, p);
            diffs_zz = zeros(p, 1);
            for k = 1:p
                [~, ~, gradk] = this.J(u, z + h(k) * v);
                fd_hv(:, k) = (gradk - grad_z) / h(k);
                diffs_zz(k) = norm(fd_hv(:, k) - Hv);
            end
            if norm(Hv) > 0
                diffs_zz = diffs_zz / norm(Hv);
            end
            disp('zz Hessian finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_zz(k))]);
            end
            disp(' ');
        end

    end
end
