%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Adv_Diff_Gaussian_Source_Objective < Dynamic_Objective

    properties
        M
        Br
        x
        time_weights
        z_time_mesh
        beta_reg
        stashed_target
        target_time_mesh
    end

    methods (Access = public)

        function [val, grad_y] = g(this, y, t)
            target = this.Evaluate_Target(t, this.x);
            val = .5 * (y - target)' * this.M * (y - target);
            grad_y = this.M * (y - target);
        end

        function [val, grad_z] = R(this, z)
            % val = .5*this.beta_reg*z'*(kron(diag(this.time_weights(2:end)),this.Br'*this.M*this.Br)*z);
            grad_z = this.beta_reg * (kron(diag(this.time_weights(2:end)), this.Br' * this.M * this.Br)) * z;
            val = .5 * z' * grad_z;
        end

        function [y_out] = g_yy_Apply(this, y_in, ~, ~)
            y_out = this.M * y_in;
        end

        function [z_out] = R_zz_Apply(this, z_in, ~)
            z_out = this.beta_reg * (kron(diag(this.time_weights(2:end)), this.Br' * this.M * this.Br)) * z_in;
        end

    end

    methods (Access = public)

        function this = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, num_space_control_nodes)
            this = this@Dynamic_Objective(n_y, n_z, T, n_t);

            % Spatial domain
            this.x = linspace(0, 1, n_y)';
            h = this.x(2) - this.x(1);

            % Mass matrix
            M = diag(4 * ones(1, n_y)) + diag(ones(1, n_y - 1), 1) + diag(ones(1, n_y - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;

            this.M = M;
            % this.beta_reg = 10^-3;
            this.beta_reg = 10^-12;

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

            load("OpInf_Training_Data.mat", "target", "target_time_mesh");
            this.stashed_target = target;
            this.target_time_mesh = target_time_mesh;
        end

        function [target] = Evaluate_Target(this, t, x)
            % target = 0.2 * t^2 * exp(-10 * (x - .5).^2);
            if t == 0
                target = zeros(size(x, 1), 1);
            else
                [~, idx] = min(abs(this.target_time_mesh - t));
                target = this.stashed_target(:, idx);
            end
        end

        function [w] = Temporal_Weights(this, t)
            w = (this.z_time_mesh - t) / (this.z_time_mesh(2) - this.z_time_mesh(1));
            Im = intersect(find(w <= 0), find(abs(w) <= 1));
            Ip = intersect(find(w > 0), find(abs(w) <= 1));
            I = find(abs(w) > 1);
            w(I) = 0;
            w(Im) = 1 + w(Im);
            w(Ip) = 1 - w(Ip);
        end

    end

end
