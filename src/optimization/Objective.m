%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Objective < handle
    % Define a scalar-valued objective function J(u,z)

    methods (Access = public)

        function this = Objective()
        end

    end

    %% Required abstract methods.

    methods (Abstract, Access = public)

        % val=J(u,z), grad_u and grad_z denote its gradients w.r.t u and z
        [val, grad_u, grad_z] = J(this, u, z)

        % Subscripts J_** denote second derivative matrices

        [u_out] = J_uu_Apply(this, u_in, u, z)

        [u_out] = J_uz_Apply(this, z_in, u, z)

        [z_out] = J_zu_Apply(this, u_in, u, z)

        [z_out] = J_zz_Apply(this, z_in, u, z)

    end

    %% Finite difference checks.

    methods (Access = public)

        function [diffs_u, diffs_z] = Finite_Difference_Gradient_Check(this, u, z)

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
