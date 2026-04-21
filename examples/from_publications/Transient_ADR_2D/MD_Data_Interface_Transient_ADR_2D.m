%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Data_Interface_Transient_ADR_2D < MD_Data_Interface

    properties
        solver
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_tmp = load('OptimizationSolution.mat', 'Y_rom').Y_rom;
            u_opt = u_tmp(:);
        end

        function [z_opt] = Load_Optimal_z(this)
            Q_rom = load('OptimizationSolution.mat', 'Q_rom').Q_rom;
            z_opt = Q_rom(:);
        end

        function [Z] = Load_Z_Data(this)
            Z = this.Load_Optimal_z();
        end

        function [D] = Load_d_Data(this)
            u_lofi = this.Load_Optimal_u();
            Y_hifi = load('OptimizationSolution.mat', 'Y_hifi').Y_hifi;
            m = size(Y_hifi, 1);
            Y_tmp = zeros(2 * m, size(Y_hifi, 3));
            Y_tmp(1:m, :) = Y_hifi(:, 1, :);
            Y_tmp((m + 1):end, :) = Y_hifi(:, 2, :);
            u_hifi = Y_tmp(:);

            D = u_hifi - u_lofi;
        end

        % Defaults to return all state elements
        % Overload function to extract component i from the state
        % Returns a vector of integers I index elements of component i
        function [I] = Separate_State_Components(this, i)
            n = this.solver.n_y / 2;
            n_t = length(this.u_opt) / (2 * n);
            if i == 1
                I = zeros(n * n_t, 1);
                for k = 1:n_t
                    I(((k - 1) * n + 1):(k * n)) = (1:n) + 2 * (k - 1) * n;
                end
            else
                I = zeros(n * n_t, 1);
                for k = 1:n_t
                    I(((k - 1) * n + 1):(k * n)) = ((n + 1):(2 * n)) + 2 * (k - 1) * n;
                end
            end
        end

        function this = MD_Data_Interface_Transient_ADR_2D(solver)
            this.solver = solver;
        end

    end

end
