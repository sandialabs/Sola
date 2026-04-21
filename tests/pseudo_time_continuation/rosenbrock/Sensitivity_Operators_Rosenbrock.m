%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Sensitivity_Operators_Rosenbrock < Euclidean_Sensitivity_Operators

    properties
        rosenbrock
    end

    methods (Access = public)

        % Reduced space objective with adjoint-based derivative
        % calculations
        function [grad, val] = Euclidean_Gradient(this, z, theta)
            [val, grad] = this.rosenbrock.J(z, theta);
        end

        function [z_out] = Euclidean_Apply_Hessian(this, z_in, z, theta)
            [~, ~, hess] = this.rosenbrock.J(z, theta);
            z_out = hess * z_in;
        end

        function [Bv] = Euclidean_Apply_B(this, theta_in, z, theta)
            B = this.rosenbrock.Compute_B(z, theta);
            Bv = B * theta_in;
        end

    end

    methods

        function this = Sensitivity_Operators_Rosenbrock(rosenbrock)
            this.rosenbrock = rosenbrock;
        end

    end
end
