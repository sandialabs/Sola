%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Dynamic_Objective < Objective
    % Define an objective function J(u,z) = \int_{0}^{T} g(y(t), t) dt + R(z),
    % represented by applying the trapezoidal rule to estimate the integral.
    % Here, the state u is the concatentation of y at all time steps.

    properties
        n_y         % Dimension of y(t).
        n_z         % Dimension of z.
        n_t         % Number of nodes in the time mesh.
        t_mesh      % Time mesh.
        w           % Quadrature weights.
    end

    properties (Dependent)
        T           % Final time.
        n_u         % Dimension n_u = n_y n_t of the full optimization state.
    end

    methods

        %% Getters for dependent properties.

        function final_time = get.T(this)
            final_time = this.t_mesh(end);
        end

        function n = get.n_u(this)
            n = this.n_y * this.n_t;
        end

        %% Constructor.

        function this = Dynamic_Objective(n_y, n_z, T, n_t)
            arguments
                n_y int32
                n_z int32
                T double
                n_t int32
            end
            this.n_y = n_y;                     
            this.n_z = n_z;                        
            this.n_t = n_t;                        
            this.t_mesh = linspace(0, T, n_t)';     
            w = ones(n_t, 1);
            w(2:end - 1) = 2;
            this.w = T * w / sum(w);                
        end

    end

    %% Required abstract methods.

    methods (Abstract, Access = public)

        [val, grad_y] = g(this, y, t)

        [val, grad_z] = R(this, z)

        [y_out] = g_yy_Apply(this, y_in, y, t)

        [z_out] = R_zz_Apply(this, z_in, z)

    end

    %% Implementation of parent class abstract methods.

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)

            val = 0;
            grad_u = 0 * u;
            for k = 1:this.n_t
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);        % y_{k} = u(I)
                [valk, gradk] = this.g(u(I), this.t_mesh(k));
                val = val + this.w(k) * valk;
                grad_u(I) = this.w(k) * gradk;
            end
            [valk, grad_z] = this.R(z);
            val = val + valk;
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            Mv = zeros(this.n_y * this.n_t, size(v, 2));
            for k = 1:this.n_t
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);        % y_{k} = u(I)
                Mv(I, :) = this.w(k) * this.g_yy_Apply(v(I, :), u(I), this.t_mesh(k));
            end
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            Mv = zeros(this.n_y * this.n_t, size(v, 2));
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            Mv = zeros(this.n_z, size(v, 2));
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            Mv = this.R_zz_Apply(v, z);
        end

    end
end
