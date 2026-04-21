%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Example_5_Constraint < Dynamic_Constraint

    % Solve the optimiation problem
    % min_{z} J(S(z),z) = int_0^T g(S(z)(t))dt + R(z)
    % where S(z) solves the ordinary differential equation
    % dy/dt = [2*t*y_1 ; 3*t^2*y_2 + t^2 - z(t)]
    % y(0) = [1, 1]
    % g(y) = (y_1-exp(t^2))^2 + (y_2-exp(t^3))^2
    % R(z) = int_0^T (z(t) - t^2)^2dt

    properties
        z_time_mesh
        weights
        beta_reg
    end

    methods (Access = public)

        function [f, f_y, f_z] = f(this, y, z, t)
            w = this.Temporal_Weights(t);
            zt = w' * z;
            f = [2 * t * y(1); 3 * t^2 * y(2) + t^2 - zt];
            f_y = [2 * t, 0; 0, 3 * t^2];
            f_z = zeros(2, this.n_z);
            f_z(2, :) = -w';
        end

        function [h, h_z] = h(this, z)
            h = [1; 1];
            h_z = zeros(2, length(z));
        end

        function [Mv] = f_yy_Apply(this, v, y, z, t, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_y, num_vecs);
        end

        function [Mv] = f_yz_Apply(this, v, y, z, t, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_y, num_vecs);
        end

        function [Mv] = f_zy_Apply(this, v, y, z, t, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_z, num_vecs);
        end

        function [Mv] = f_zz_Apply(this, v, y, z, t, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_z, num_vecs);
        end

        function [Mv] = h_zz_Apply(this, v, z, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_z, num_vecs);
        end

    end

    methods (Access = public)

        function this = Example_5_Constraint(n_y, n_z, T, n_t)
            this = this@Dynamic_Constraint(n_y, n_z, T, n_t);
            this.z_time_mesh = linspace(0, T, n_z + 1)';
            this.z_time_mesh = this.z_time_mesh(2:end);

            weights = ones(this.n_z + 1, 1);
            weights(1) = .5;
            weights(end) = .5;
            weights = this.t_mesh(end) * weights / sum(weights);
            weights = weights(2:end);
            this.weights = weights;
            this.beta_reg = 10^-4;
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
