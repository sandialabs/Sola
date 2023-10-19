% Define the objective function J(u, z) where
% $u \in R^{n_u}$
% $z \in R^{n_z}$
% $J(u, z) \in R$

classdef Objective < handle

    methods (Abstract, Access = public)

         % Input:
         % u: the state u in R^{n_u}
         % z: the control z \in R^{n_z}
         % Output:
         % val: J(u, z) in R
         % grad_u: \nabla_u J(u, z) in R^{n_u}
         % grad_z: \nabla_z J(u, z) in R^{n_z}
        [val, grad_u, grad_z] = J(this, u, z);

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % Mv: \nabla_{u, u} J(u, z)v in R^{n_u}
        [Mv] = J_uu_Apply(this, v, u, z);

        % Input:
        % v: a direction v in R^{n_z}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % Mv: \nabla_{u, z} J(u, z)v in R^{n_u}
        [Mv] = J_uz_Apply(this, v, u, z);

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output
        % Mv: \nabla_{z, u} J(u, z)v in R^{n_z}
        [Mv] = J_zu_Apply(this, v, u, z);

        % Input:
        % v: a direction v in R^{n_z}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output
        % Mv: \nabla_{z, z} J(u, z)v in R^{n_z}
        [Mv] = J_zz_Apply(this, v, u, z);

    end

    methods (Access = public)

        function this = Objective( )

        end

        % Finite difference test functions

        % Input:
        % u: the control u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % diffs_u: vector of finite difference errors for u gradient
        % diffs_z: vector of finite difference errors for z gradient
        function [diffs_u, diffs_z] = Finite_Difference_Gradient_Check(this, u, z)
            [val, grad_u, grad_z] = this.J(u, z);

            h = 10.^(-2:-1:-6);
            p = length(h);

            m = length(grad_u);
            du = randn(m, 1);
            du = du / norm(du);
            grad_du = du'*grad_u;
            fd_grad = zeros(p, 1);
            for k = 1:p
                valk = this.J(u + h(k)*du, z);
                fd_grad(k) = (valk - val)/h(k);
            end

            diffs_u = abs(grad_du - fd_grad) / abs(grad_du);
            disp('u gradient finite difference check')
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_u(k))])
            end
            disp(' ')

            n = length(grad_z);
            dz = randn(n, 1);
            dz = dz / norm(dz);
            grad_dz = dz' * grad_z;
            fd_grad = zeros(p, 1);
            for k = 1:p
                valk = this.J(u, z + h(k)*dz);
                fd_grad(k) = (valk - val) / h(k);
            end

            diffs_z = abs(grad_dz - fd_grad) / abs(grad_dz);
            disp('z gradient finite difference check')
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_z(k))])
            end
            disp(' ')


        end

        % Input:
        % u: the control u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % diffs_uu: vector of finite difference errors for uu hessian
        % diffs_uz: vector of finite difference errors for uz hessian
        % diffs_zu: vector of finite difference errors for zu hessian
        % diffs_zz: vector of finite difference errors for zz hessian
        function [diffs_uu, diffs_uz, diffs_zu, diffs_zz] = Finite_Difference_Hessian_Check(this, u, z)
            [~, grad_u, grad_z] = this.J(u, z);

            h = 10.^(-2:-1:-6);
            p = length(h);
            m = length(grad_u);
            n = length(grad_z);

            v = randn(m, 1);
            v = v/norm(v);
            Hv = this.J_uu_Apply(v, u, z);
            fd_hv = zeros(m, p);
            diffs_uu = zeros(p, 1);
            for k = 1:p
                [~, gradk] = this.J(u + h(k)*v, z);
                fd_hv(:, k) = (gradk - grad_u) / h(k);
                diffs_uu(k) = norm(fd_hv(:, k) - Hv);
            end
            if norm(Hv) > 0
                diffs_uu = diffs_uu / norm(Hv);
            end
            disp('uu Hessian finite difference check')
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_uu(k))])
            end
            disp(' ')

            v = randn(n, 1);
            v = v/norm(v);
            Hv = this.J_uz_Apply(v, u, z);
            fd_hv = zeros(m, p);
            diffs_uz = zeros(p, 1);
            for k = 1:p
                [~, gradk] = this.J(u, z + h(k)*v);
                fd_hv(:, k) = (gradk - grad_u) / h(k);
                diffs_uz(k) = norm(fd_hv(:, k) - Hv);
            end
            if norm(Hv) > 0
                diffs_uz = diffs_uz / norm(Hv);
            end
            disp('uz Hessian finite difference check')
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_uz(k))])
            end
            disp(' ')

            v = randn(m, 1);
            v = v/norm(v);
            Hv = this.J_zu_Apply(v, u, z);
            fd_hv = zeros(n, p);
            diffs_zu = zeros(p, 1);
            for k = 1:p
                [~, ~, gradk] = this.J(u + h(k)*v, z);
                fd_hv(:, k) = (gradk - grad_z) / h(k);
                diffs_zu(k) = norm(fd_hv(:, k) - Hv);
            end
            if norm(Hv) > 0
                diffs_zu = diffs_zu / norm(Hv);
            end
            disp('zu Hessian finite difference check')
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_zu(k))])
            end
            disp(' ')

            v = randn(n, 1);
            v = v / norm(v);
            Hv = this.J_zz_Apply(v, u, z);
            fd_hv = zeros(n, p);
            diffs_zz = zeros(p, 1);
            for k = 1:p
                [~, ~, gradk] = this.J(u, z + h(k)*v);
                fd_hv(:, k) = (gradk - grad_z) / h(k);
                diffs_zz(k) = norm(fd_hv(:, k) - Hv);
            end
            if norm(Hv)>0
                diffs_zz = diffs_zz / norm(Hv);
            end
            disp('zz Hessian finite difference check')
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_zz(k))])
            end
            disp(' ')
        end

    end
end
