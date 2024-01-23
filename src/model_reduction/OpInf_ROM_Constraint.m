classdef OpInf_ROM_Constraint < Dynamic_Constraint
    % Operator Inference Reduced-order Model.
    %
    % This class represents systems of ordinary differential equations (ODEs)
    % with the polynomial
    %
    % .. math:: \ddt\y(t) = \c + \A\y(t) + \H[\y(t)\otimes\y(t)] + \B\q(t),
    %
    % where
    %
    % * :math:`\y(t)\in\R^{n_y}` is the ODE state,
    % * :math:`\q(t)\in\R^{n_q}` is the time-dependent (state-independent) input, and
    % * :math:`\c,\A,\H,\B` are the system "operators".
    %
    % The system operators are determined by solving a regression problem,
    %
    % .. math::
    %  \min_{\c,\A,\H,\B}\sum_{j=1}^{n_t}\left\|
    %  \c + \A\y_j + \H[\y_j\otimes\y_j] + \B\q_j - \dot{\y}_j
    %  \right\|_{2}^{2}
    %  + \left\|\bfGamma[~\c~~\A~~\H~~\B~]\trp\right\|_{F}^{2},
    %
    % where
    %
    % * :math:`\y_j\in\R^{n_y}` is a snapshot :math:`\y(t_j)` of the ODE state at time :math:`t_j`,
    % * :math:`\q_j\in\R^{n_q}` is the input :math:`\q(t_j)` at time :math:`t_j`,
    % * :math:`\dot{\y}_j\in\R^{n_y}` is the time derivative of the ODE state at time :math:`t_j`,
    %   i.e., :math:`\ddt\y(t)\big|_{t = t_j}`,
    % * :math:`n_t \in \N` is the number of snapshots, and
    % * :math:`\bfGamma\in\R^{d \times d}` is a Tikhonov regularization matrix.
    %
    % This class inherits from :class:`Dynamic_Constraint`, hence the model is integrated in time
    % using the implicit Euler method. In this context, the optimization state and control are
    % given by
    %
    % .. math::
    %  \u = \left(\begin{array}{c}\y_1 \\ \vdots \\ \y_{n_t}\end{array}\right)\in\R^{n_u},
    %  \qquad
    %  \z = \left(\begin{array}{c}\q_1 \\ \vdots \\ \q_{n_t}\end{array}\right)\in\R^{n_z},
    %
    % with :math:`n_u = n_y n_t` and :math:`n_z = n_q n_t`.
    % Note that the state and the control share a time mesh.
    %
    % The user indicates the terms to appear in the system.
    %
    % Example
    % -------
    % .. code-block:: matlab
    %
    %  y0 = randn(20, 1);
    %  operators = {Constant_Operator(), Linear_Operator(), ...
    %               Quadratic_Operator(), Input_Operator()};
    %  constraint = OpInf_ROM_Constraint(20, 4, 1, 400, y0, operators);

    properties
        regularizer         % Tikhonov regularization matrix.
        operators           % Operator objects defining the differential system.
    end

    properties (SetAccess = protected)
        n_q                 % Dimension :math:`n_q` of the input :math:`\q(t)`.
        ds
        d                   % Number of operator entries for each system mode.
        y0                  % Initial condition :math:`\y(0)\in\R^{n_y}` for the ODE.
    end

    methods

        %% Property setters.

        function set.operators(this, operators)
            % Set the dynamical system operators. These may be uninitialized,
            % meaning Set_Operators() has not been called.
            num_operators = length(operators);
            this.ds = zeros(num_operators, 1);
            for i = 1:num_operators
                this.ds(i) = operators{i}.Column_Dimension(this.n_y, this.n_q);
            end
            this.d = sum(this.ds);
            this.operators = operators;
        end

        function set.regularizer(this, reg)
            % Set the Tikhonov regularizer.
            %
            % Parameters
            % ----------
            % reg
            %   Tikhonov regularization, one of the following:
            %
            %   * scalar: regularizer = reg * I
            %   * vector of length ``length(this.operators)``:
            %     regularizer = diag(reg(1)I_1, reg(2)I_2, ...) where each I_i
            %     is the column dimension of the ith operator.
            %   * vector of length d: regularizer = diag(reg)
            %   * matrix of size (d, d): regularizer = reg.

            % Scalar regularizer.
            if isscalar(reg)
                % Scalar regularizer --> Gamma = regularizer * I
                reg = reg * eye(this.d);
            end

            % Vector regularizer.
            num_operators = length(this.operators);
            if min(size(reg)) == 1
                % Vector regularizer
                if length(reg) == num_operators
                    % Scalar regularizer for each operator
                    newreg = [];
                    for i = 1:num_operators
                        newreg = horzcat(newreg, reg(i) * ones(1, this.ds(i)));
                    end
                    reg = newreg;
                end
                reg = diag(reg);
            end

            % Matrix regularizer.
            if size(reg, 1) ~= this.d || size(reg, 2) ~= this.d
                error('invalid regularizer');
            end

            this.regularizer = reg;
        end

    end

    methods (Access = public)

        %% Constructor.

        function this = OpInf_ROM_Constraint(n_y, n_q, T, n_t, y0, operators)
            % Parameters
            % ----------
            % n_y
            %   Dimension :math:`n_y` of the state :math:`\y_{j}` at each time.
            %   Note that if this is a reduced-order model, this is the REDUCED
            %   state dimension, denoted elsewhere as :math:`n_y'`.
            % n_q
            %   Dimension :math:`n_q` of the input :math:`\q`.
            % T
            %   Final time :math:`T > 0`.
            % n_t
            %   Number of nodes :math:`N` in the time mesh.
            % y0
            %   Initial condition :math:`\y(0)\in\R^{n_y}` for the ODE.
            % operators
            %   ODE system operators. These may be uninitialized,
            %   i.e., with null entries.

            this@Dynamic_Constraint(n_y, n_q * n_t, T, n_t);
            this.n_q = n_q;
            this.y0 = y0;
            this.operators = operators;
            this.regularizer = 0;
        end

        %% Training.

        function Learn_Operators(this, Y, Q, dYdt)
            % Infer the entries of the system operators by solving the
            % regression problem
            %
            % .. math::
            %  \min_{\c,\A,\H,\B}\sum_{j=1}^{n_t}\left\|
            %  \c + \A\y_j + \H[\y_j\otimes\y_j] + \B\q_j - \dot{\y}_j
            %  \right\|_{2}^{2}
            %  + \left\|\bfGamma[~\c~~\A~~\H~~\B~]\trp\right\|_{F}^{2}
            %
            % where
            %
            % * :math:`\Y = [~\y_1~~\cdots~~\y_{n_t}~] \in \R^{n_y \times n_t}`
            %   are state snapshots,
            % * :math:`\dot{\Y} = [~\dot{\y}_1~~\cdots~~\dot{\y}_{n_t}~]
            %   \in \R^{n_y \times n_t}` collects the time derivative of the
            %   state snapshots,
            % * :math:`\Q = [~\q_1~~\cdots~~\q_{n_t}~] \in \R^{n_q \times n_t}`
            %   are the corresponding inputs, and
            % * :math:`\bfGamma\in\R^{d \times d}` is a Tikhonov regularization
            %   matrix (the ``regularizer`` property).
            %
            % Parameters
            % ----------
            % Y
            %   State snapshots :math:`\Y \in \R^{n_y \times n_t}`.
            % Q
            %   Inputs :math:`\Q \in \R^{n_q \times n_t}`.
            % dYdt
            %   (Optional) Time derivative of snapshots
            %   :math:`\dot{\Y} \in \R^{n_y \times n_t}`.
            %   If not provided, these are estimated from ``Y`` with finite differences.
            arguments
                this
                Y (:, :) {mustBeNumeric}
                Q (:, :) {mustBeNumeric}
                dYdt (:, :) {mustBeNumeric} = []
            end

            % If time derivatives not given, estimate with backward Euler (to match time stepping scheme).
            if size(dYdt, 1) == 0
                dt = this.t_mesh(2) - this.t_mesh(1);
                dYdt = (Y(:, 2:end) - Y(:, 1:(end - 1))) / dt;
                Y = Y(:, 2:end);
                Q = Q(:, 2:end);
            end

            % Construct the data matrix.
            D = [];
            for i = 1:length(this.operators)
                D = horzcat(D, this.operators{i}.Datablock(Y, Q)');
            end

            % disp(['Data matrix size: (', num2str(size(D, 1)), ', ', num2str(size(D, 2)), ')']);
            % disp(['RHS matrix size: (', num2str(size(dYdt, 2)), ', ', num2str(size(dYdt, 1)), ')']);

            % Solve the regression problem for the operator entries.
            if norm(this.regularizer) > 0
                % Tikhonov regularization.
                Ohat = lscov([D; this.regularizer], [dYdt'; zeros(this.d, this.n_y)])';
            else
                % No regularization.
                Ohat = lscov(D, dYdt')';
            end

            % disp(['Operator matrix size: (', num2str(size(Ohat, 1)), ', ', num2str(size(Ohat, 2)), ')']);

            % Unpack the operator entries.
            index = 1;
            for i = 1:length(this.operators)
                newindex = index + this.ds(i);
                this.operators{i}.Set_Entries(Ohat(:, index:(newindex - 1)));
                index = newindex;
            end
        end

        %% Implement abstract methods from the parent class.

        function [f, f_y, f_z] = f(this, y, z, t)

            % Extract the input q(t) from the control z based on the time t.
            I = this.Input_Indices(t);
            q = z(I);

            % Allocate space for outputs.
            f = zeros(this.n_y, 1);
            f_y = zeros(this.n_y, this.n_y);
            f_z = zeros(this.n_y, this.n_z);

            % Evaluate each operator.
            for i = 1:length(this.operators)
                op = this.operators{i};
                f = f + op.Apply(y, q);
                f_y = f_y + op.Jacobian_y(y, q);
                f_z(:, I) = f_z(:, I) + op.Jacobian_q(y, q);
            end
        end

        function [h, h_z] = h(this, ~)
            h = this.y0;
            h_z = zeros(this.n_z, 1);
        end

        function [Mv] = f_yy_Apply(this, v, y, z, t, lambda)
            I = this.Input_Indices(t);
            q = z(I);
            Mv = zeros(this.n_y, size(v, 2));
            for i = 1:length(this.operators)
                Mv = Mv + this.operators{i}.Hessian_yy_Apply(v, y, q, lambda);
            end
        end

        function [Mv] = f_yz_Apply(this, v, y, z, t, lambda)
            I = this.Input_Indices(t);
            q = z(I);
            vt = v(I, :);
            Mv = zeros(this.n_y, size(v, 2));
            for i = 1:length(this.operators)
                Mv = Mv + this.operators{i}.Hessian_yq_Apply(vt, y, q, lambda);
            end
        end

        function [Mv] = f_zy_Apply(this, v, y, z, t, lambda)
            I = this.Input_Indices(t);
            q = z(I);
            Mv = zeros(this.n_z, size(v, 2));
            for i = 1:length(this.operators)
                Mv(I, :) = Mv(I, :) + this.operators{i}.Hessian_qy_Apply(v, y, q, lambda);
            end
        end

        function [Mv] = f_zz_Apply(this, v, y, z, t, lambda)
            I = this.Input_Indices(t);
            q = z(I);
            vt = v(I, :);
            Mv = zeros(this.n_z, size(v, 2));
            for i = 1:length(this.operators)
                Mv(I, :) = Mv(I, :) + this.operators{i}.Hessian_yq_Apply(vt, y, q, lambda);
            end
        end

        function [Mv] = h_zz_Apply(this, v, ~, ~)
            Mv = zeros(this.n_z, size(v, 2));
        end

    end

    methods (Access = protected)

        function [mask] = Input_Indices(this, t)
            [~, t_index] = min(abs(t - this.t_mesh));
            mask = (this.n_q * (t_index - 1) + 1):(this.n_q * t_index);
        end

    end
end
