classdef MD_Laplacian_Like_Operator_Properties < handle

    properties

    end

    methods

        function this = MD_Laplacian_Like_Operator_Properties()

        end

        function [val] = Get_Rectangular_Domain_Squared_Inv_Operator_Trace(this, beta, nodes)

            if size(nodes, 2) == 1
                Lx = max(nodes(:, 1)) - min(nodes(:, 1));
                n = length(nodes(:, 1)) - 1;
                e = 1 + beta * (pi / Lx)^2 * (0:n).^2;
                e = e';
            elseif size(nodes, 2) == 2
                Lx = max(nodes(:, 1)) - min(nodes(:, 1));
                Ly = max(nodes(:, 2)) - min(nodes(:, 2));
                n = round(sqrt(length(nodes(:, 1)))) - 1;
                e = 1 + beta * pi^2 * (kron(((0:n).^2)', ones(n + 1, 1)) / Lx^2 + kron(ones(n + 1, 1), ((0:n).^2)') / Ly^2);
            else
                disp('Get_Squared_Inv_Operator_Trace error: Dimensions greater than 2 are not supported.');
            end

            val = sum((1 ./ e).^2);
        end

        function [val] = Randomized_Inv_Operator_Trace_Estimation(this, u_prior_interface, num_samples)
            u_out = u_prior_interface.Sample_with_Covariance_W_u_Acute_Inverse(num_samples);
            tmp = u_prior_interface.Apply_M_u(u_out);
            val = (1 / num_samples) * sum(diag(u_out' * tmp));
        end

    end
end
