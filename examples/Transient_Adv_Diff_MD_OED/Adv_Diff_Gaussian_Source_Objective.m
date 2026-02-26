classdef Adv_Diff_Gaussian_Source_Objective < Dynamic_Objective

    properties
        M
        Br
        x
        time_weights
        z_time_mesh
        reg_coeff
    end

    methods (Access = public)

        function [val, grad_y] = g(this, y, t)
            target = this.Evaluate_Target(t, this.x);
            val = .5 * (y - target)' * this.M * (y - target);
            grad_y = this.M * (y - target);
        end

        function [val, grad_z] = R(this, z)
            grad_z = this.reg_coeff * (kron(diag(this.time_weights(2:end)), this.Br' * this.M * this.Br)) * z;
            val = .5 * z' * grad_z;
        end

        function [Mv] = g_yy_Apply(this, v, ~, ~)
            Mv = this.M * v;
        end

        function [Mv] = R_zz_Apply(this, v, ~)
            Mv = this.reg_coeff * (kron(diag(this.time_weights(2:end)), this.Br' * this.M * this.Br)) * v;
        end

    end

    methods (Access = public)

        function this = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, num_space_control_nodes, reg_coeff)
            this = this@Dynamic_Objective(n_y, n_z, T, n_t);

            % Spatial domain
            this.x = linspace(0, 1, n_y)';
            h = this.x(2) - this.x(1);

            M = diag(4 * ones(1, n_y)) + diag(ones(1, n_y - 1), 1) + diag(ones(1, n_y - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;

            this.M = M;
            this.reg_coeff = reg_coeff;

            % Trapezoid rule for time integration
            time_weights = ones(n_t, 1);
            time_weights(1) = .5 * time_weights(1);
            time_weights(end) = .5 * time_weights(end);
            time_weights = T * time_weights / sum(time_weights);
            this.time_weights = time_weights;

            % Control matrix (Gaussian)
            source_loc = linspace(0, 1, num_space_control_nodes)';
            Br = zeros(n_y, num_space_control_nodes);
            for k = 1:num_space_control_nodes
                Br(:, k) = exp(-200 * (this.x - source_loc(k)).^2);
            end
            this.Br = Br;
            this.z_time_mesh = linspace(0, T, n_t)';
            this.z_time_mesh = this.z_time_mesh(2:end);
        end

        function [target] = Evaluate_Target(this, t, x)
            target = 0.2 * t^2 * exp(-10 * (x - .5).^2);
        end

    end

end
