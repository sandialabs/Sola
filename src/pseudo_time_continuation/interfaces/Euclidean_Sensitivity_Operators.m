%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Euclidean_Sensitivity_Operators < Sensitivity_Operators

    properties

    end

    methods (Abstract, Access = public)

        [grad, val] = Euclidean_Gradient(this, z, theta)

        [z_out] = Euclidean_Apply_Hessian(this, z_in, z, theta)

        [z_out] = Euclidean_Apply_B(this, theta_in, z, theta)

    end

    methods

        function this = Euclidean_Sensitivity_Operators()

        end

        function [grad, val] = Gradient(this, z, theta_traj, time_index)
            [grad, val] = this.Euclidean_Gradient(z, theta_traj.Get_theta_n(time_index));
        end

        function [z_out] = Apply_Hessian(this, z_in, z, theta_traj, time_index)
            z_out = this.Euclidean_Apply_Hessian(z_in, z, theta_traj.Get_theta_n(time_index));
        end

        function [z_out] = Apply_B(this, z, theta_traj, time_index)
            z_out = this.Euclidean_Apply_B(theta_traj.Get_dtheta(), z, theta_traj.Get_theta_n(time_index));
        end

        function [] = Finite_Difference_Gradient_Check(this, z, theta)
            [grad, val] = this.Euclidean_Gradient(z, theta);
            h = 10.^(-1:-1:-6);
            v = randn(length(z), 1);
            error = zeros(length(h), 1);
            for k = 1:length(h)
                [~, valk] = this.Euclidean_Gradient(z + h(k) * v, theta);
                error(k) = abs((valk - val) / h(k) - grad' * v) / abs(grad' * v);
            end

            disp('Gradient finite difference test');
            for k = 1:length(h)
                disp(['Step size = ', num2str(h(k)), ' and error = ', num2str(error(k))]);
            end
        end

        function [diffs] = Finite_Difference_Hessian_Check(this, z, theta)
            grad = this.Euclidean_Gradient(z, theta);
            n = length(grad);
            v = randn(n, 1);
            v = v / norm(v);
            Hv = this.Euclidean_Apply_Hessian(v, z, theta);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_hv = zeros(n, p);
            diffs = zeros(p, 1);
            for k = 1:p
                gradk = this.Euclidean_Gradient(z + h(k) * v, theta);
                fd_hv(:, k) = (gradk - grad) / h(k);
                diffs(k) = norm(fd_hv(:, k) - Hv) / norm(Hv);
            end
            disp('Hessian finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
            end
            disp(' ');
        end

        function [diffs] = Finite_Difference_B_Check(this, z, theta)
            grad = this.Euclidean_Gradient(z, theta);
            n = length(theta);
            v = randn(n, 1);
            v = v / norm(v);
            Bv = this.Euclidean_Apply_B(v, z, theta);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_Bv = zeros(length(z), p);
            diffs = zeros(p, 1);
            for k = 1:p
                gradk = this.Euclidean_Gradient(z, theta + h(k) * v);
                fd_Bv(:, k) = (gradk - grad) / h(k);
                diffs(k) = norm(fd_Bv(:, k) - Bv) / norm(Bv);
            end
            disp('B finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
            end
            disp(' ');
        end

    end
end
