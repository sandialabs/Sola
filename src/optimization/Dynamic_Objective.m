classdef Dynamic_Objective < Objective
    % Define an objective function
    %
    % .. math:: J(\u,\z) = \int_{0}^{T} g(\u(t), t) dt + R(\z),
    %
    % represented by applying the trapezoidal rule to estimate the integral:
    %
    % .. math:: J(\u,\z) = \sum_{j=1}^{N} w_{j} g(\y_{j}, t_{j}) + R(\z).
    %
    % Here,
    %
    % * :math:`\u = (\y_{1}\trp \cdots \y_{n_t}\trp)\trp\in\R^{n_u}`
    %   with :math:`\y_{j} \in \R^{n_y}` and :math:`n_u = n_y n_t`,
    % * :math:`0 = t_{1} < t_{2} < \cdots < t_{n_t} = T`
    %   with equal spacing :math:`\delta t = t_{j+1} - t_{j}`,
    % * :math:`\z \in \R^{n_z}`,
    % * :math:`g:\R^{n_y} \times [0, T] \to \R`,
    % * :math:`R:\R^{n_z} \to \R`, and
    % * :math:`w_{j} = \delta t` for :math:`j = 2, \ldots, N - 1`
    %   and :math:`w_{1} = w_{n_t} = \frac{1}{2}\delta t`.

    properties
        n_y         % Dimension :math:`n_y` of the differential equation state :math:`\y`.
        n_z         % Dimension :math:`n_z` of the control :math:`\z`.
        n_t         % Number of nodes :math:`n_t` in the time mesh.
        t_mesh      % Time mesh :math:`(t_1,\ldots,t_{n_t})\trp`.
        w           % Quadrature weights :math:`(w_1,\ldots,w_{n_t})\trp` for the time integral.
    end

    properties (Dependent)
        T           % Final time :math:`T = t_{n_t}`.
    end

    methods

        function final_time = get.T(this)
            final_time = this.t_mesh(end);
        end

    end

    methods (Abstract, Access = public)

        [val, grad_y] = Time_Instance_Objective(this, y, t)
        % Evaluate the integrand :math:`g(\y,t)`
        % and its gradient :math:`\grad{y}g(\y,t)`.
        %
        % Parameters
        % ----------
        % y
        %   Differential equation state :math:`\y\in\R^{n_y}`.
        % t
        %   Time :math:`t`.
        %
        % Returns
        % -------
        % val : double
        %   Function value :math:`g(\y,t)\in\R`.
        % grad_y : vector
        %   Function gradient :math:`\grad{y}g(\y,t)\in\R^{n_y}`.

        [val, grad_z] = Regularization_Objective(this, z)
        % Evaluate the regularization term :math:`R(\z)`
        % and its gradient :math:`\grad{z}R(\z)`.
        %
        % Parameters
        % ----------
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        %
        % Returns
        % -------
        % val : double
        %   Function value :math:`R(\z)\in\R`.
        % grad_y : vector
        %   Function gradient :math:`\grad{z}R(\z)\in\R^{n_z}`.

        [Mv] = Time_Instance_Objective_yy_Apply(this, v, y, t)
        % Compute the Hessian-vector product :math:`\grad{y,y}g(\y, t)\v`.
        %
        % Parameters
        % ----------
        % v
        %   Search direction :math:`\v\in\R^{n_y}`.
        % y
        %   Differential equation state :math:`\y\in\R^{n_y}`.
        % t
        %   Time :math:`t`.
        %
        % Returns
        % -------
        % Mv : vector
        %   Hessian-vector product :math:`\grad{y,y}g(\y, t)\v\in\R^{n_y}`

        [Mv] = Regularization_Objective_zz_Apply(this, v, z)
        % Compute the Hessian-vector product :math:`\grad{z,z}R(\z)\v`.
        %
        % Parameters
        % ----------
        % v
        %   Search direction :math:`\v\in\R^{n_z}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        %
        % Returns
        % -------
        % Mv : vector
        %   Hessian-vector product :math:`\grad{z,z}R(\z)\v\in\R^{n_z}`

    end

    methods (Access = public)

        %% Implementation of parent class abstract methods.

        function [val, grad_u, grad_z] = J(this, u, z)
            % Evaluate the objective function and its derivatives.
            %
            % Parameters
            % ----------
            % u
            %   State :math:`\u\in\R^{n_u}`.
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            %
            % Returns
            % -------
            % val : double
            %   Objective value :math:`J(\u,\z) \in \R`.
            % grad_u : vector
            %   Objective gradient :math:`\grad{u}J(\u,\z) \in \R^{n_u}`.
            % grad_z : vector
            %   Objective gradient :math:`\grad{z}J(\u,\z) \in \R^{n_z}`.

            val = 0;
            grad_u = 0 * u;
            for k = 1:this.n_t
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);        % y_{k} = u(I)
                [valk, gradk] = this.Time_Instance_Objective(u(I), this.t_mesh(k));
                val = val + this.w(k) * valk;
                grad_u(I) = this.w(k) * gradk;
            end
            [valk, grad_z] = this.Regularization_Objective(z);
            val = val + valk;
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            Mv = zeros(this.n_y * this.n_t, size(v, 2));
            for k = 1:this.n_t
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);        % y_{k} = u(I)
                Mv(I, :) = this.w(k) * this.Time_Instance_Objective_yy_Apply(v(I, :), u(I), this.t_mesh(k));
            end
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            Mv = zeros(this.n_y * this.n_t, size(v, 2));
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            Mv = zeros(this.n_z, size(v, 2));
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            Mv = this.Regularization_Objective_zz_Apply(v, z);
        end

        %% Constructor.

        function this = Dynamic_Objective(n_y, n_z, T, n_t)
            % Parameters
            % ----------
            % n_y : int
            %   Dimension :math:`n_y` of the differential equation state :math:`\y`.
            % n_z : int
            %   Dimension :math:`n_z` of the control :math:`\z`.
            % T : double
            %   Final time :math:`T`.
            % n_t : int
            %   Number of nodes :math:`n_t` in the time mesh.
            this.n_y = n_y;                         % ODE state dimension
            this.n_z = n_z;                         % Control dimension
            this.n_t = n_t;                         % Number of time nodes
            this.t_mesh = linspace(0, T, n_t)';     % Discrete time domain
            w = ones(n_t, 1);
            w(2:end - 1) = 2;
            this.w = T * w / sum(w);                % Quadrature weights for time integral
        end

    end

end
