classdef Tracer_Optimization < handle
    properties
        obj                 % Instance of a subclass of :class:`Objective`.
        con                 % Instance of a subclass of :class:`Constraint`.
    end

    methods (Access = public)

        function this = Tracer_Optimization(obj, con)
            this.obj = obj;
            this.con = con;
        end

        %% Optimization functions

        function [u, z] = Optimize(this, z0)

        end

        function [val, grad, hessian_data] = Jhat(this, z)
            u = this.con.State_Solve(z);
            [val, grad_u, grad_z] = this.obj.J(u, z);
            lambda = this.con.c_u_Transpose_Inverse_Apply(-grad_u, u, z);
            grad = this.con.c_z_Transpose_Apply(lambda, u, z);
            grad = grad + grad_z;
            hessian_data = [u; z; lambda];
        end

        function [z_out] = Jhat_hessVec(this, hessian_data, z_in)
            % Compute the vector-Hessian-vector product
            % :math:`\bflambda\trp\grad{z,z}\hat{J}(\z)\v`
            % via an adjoint-based approach.
            %
            % Parameters
            % ----------
            % hessian_data : vector
            %   Concatenation of the state :math:`\u\in\R^{n_u}`,
            %   control :math:`\z\in\R^{n_z}`,
            %   and adjoint :math:`\bflambda\in\R^{n_u}`.
            % z_in : vector
            %   Control direction :math:`\v\in\R^{n_z}`.
            %
            % Returns
            % -------
            % z_out : vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp\grad{z,z}\hat{J}(\z)\v\in\R^{n_z}`.

            % Extract state, control, and adjoint from hessian_data.

        end

        %% Finite difference tests

        function [diffs] = Finite_Difference_Gradient_Check(this, z)
            % Check the implementation of :meth:`Jhat()` via finite differences.
            %
            % Parameters
            % ----------
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            %
            % Returns
            % -------
            % diffs : vector
            %   Finite difference errors.

            [val, grad] = this.Jhat(z);
            n = length(grad);
            dz = randn(n, 1);
            dz = dz / norm(dz);
            grad_dz = dz' * grad;
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_grad = zeros(p, 1);
            for k = 1:p
                valk = this.Jhat(z + h(k) * dz);
                fd_grad(k) = (valk - val) / h(k);
            end

            diffs = abs(grad_dz - fd_grad) / abs(grad_dz);
            if this.verbose
                disp('Gradient finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

        function [diffs] = Finite_Difference_Hessian_Check(this, z)
            % Check the implementation of :meth:`Jhat_hessVec()` via finite differences.
            %
            % Parameters
            % ----------
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            %
            % Returns
            % -------
            % diffs : vector
            %   Finite difference errors.

            [~, grad, hessian_data] = this.Jhat(z);
            n = length(grad);
            v = randn(n, 1);
            v = v / norm(v);
            Hv = this.Jhat_hessVec(hessian_data, v);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_hv = zeros(n, p);
            diffs = zeros(p, 1);
            for k = 1:p
                [~, gradk] = this.Jhat(z + h(k) * v);
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

    end
end
