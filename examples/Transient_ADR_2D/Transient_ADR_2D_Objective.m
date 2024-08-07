classdef Transient_ADR_2D_Objective < Dynamic_Objective
    % Objective function for a two-species time-dependent
    % advection-diffusion-reaction problem in two spatial dimensions.
    % This objective measures the relative concentration of a
    % contaminant near a specified protection zone. Mathematically, the
    % objective emulates the function
    %
    % .. math::
    %    \frac{1}{2}\int_0^T \|\y_1(t) \ast \p\|_{\M}^{2}\dt
    %    + \frac{\gamma}{2}\int_{t_2}^T \|\q(t)\|_{2}^{2}\dt
    %
    % where :math:`\p \in \R^{n_x}` is a discretization of a wide Gaussian
    % function centered over a protection zone and in which :math:`\ast`
    % denotes the Hadamard (element-wise) product.

    properties
        beta_reg        % Regularization hyperparameter for the control.
    end

    properties (SetAccess = protected)
        x               % :math:`x`-coordinates of the spatial nodes.
        y               % :math:`y`-coordinates of the spatial nodes.
        n_x             % Number of spatial nodes.
        M               % Mass matrix corresponding to the spatial coordinates.
        n_q             % Number of control nodes at any fixed time.
        % Bq              % Control nodes spatial discretization.
        w_z             % Weights for the time integral in the control regularization.
        target_weight   % Weight vector :math:`\p` for measuring contaminant in target areas.
    end

    properties (Access = protected)
        pMp             % Mass matrix with target weights broadcast to rows and columns.
    end

    methods (Access = public)

        function [val, grad_y] = g(this, y, ~)
            % :math:`|| y_1(t) .* p ||_M^2 / 2`
            y1 = y(1:this.n_x);
            p = this.target_weight;
            grad_y1 = p .* (this.M * (y1 .* p));
            grad_y2 = zeros(this.n_x, 1);
            val = .5 * y1' * grad_y1;
            grad_y = [grad_y1; grad_y2];
        end

        function [val, grad_z] = R(this, z)
            % :math:`\gamma \int_{t_2}^{T} || q(t) ||_2^2 dt / 2`
            val = .5 * this.beta_reg * this.w_z' * (z.^2);
            grad_z = this.beta_reg * this.w_z .* z;
        end

        function [y_out] = g_yy_Apply(this, y_in, ~, ~)
            y_in1 = y_in(1:this.n_x, :);
            y_out = zeros(this.n_y, size(y_in, 2));
            y_out(1:this.n_x, :) = this.pMp * y_in1;
        end

        function [z_out] = R_zz_Apply(this, z_in, ~)
            z_out = this.beta_reg * this.w_z .* z_in;
        end

    end

    methods (Access = public)

        function this = Transient_ADR_2D_Objective(center, x, y, M, T, n_t, n_q, beta_reg)
            arguments
                center {mustBeNumeric}
                x {mustBeNumeric}
                y {mustBeNumeric}
                M {mustBeNumeric}
                T {mustBePositive}
                n_t {mustBePositive, mustBeInteger}
                n_q {mustBePositive, mustBeInteger}
                % Bq {mustBeNumeric}
                beta_reg {mustBePositive} = 1e-3
            end

            % Dimensions.
            n_x = numel(x);
            n_y = 2 * n_x;
            % n_q = size(Bq, 2);
            n_z = n_q * (n_t - 1);
            this = this@Dynamic_Objective(n_y, n_z, T, n_t);

            % Spatial properties.
            this.n_x = n_x;
            this.x = reshape(x, n_x, 1);
            this.y = reshape(y, n_x, 1);
            if (size(M, 1) ~= n_x) || (size(M, 2) ~= n_x)
                error('Mass matrix not aligned with spatial coordinates');
            end
            this.M = M;
            % Target weight.
            this.target_weight = exp(-10 .* sum(([this.x, this.y] - reshape(center, 1, [])).^2, 2));
            this.pMp = this.target_weight' .* this.M .* this.target_weight;

            % Control properties
            this.n_q = n_q;
            % if size(this.Bq, 1) ~= n_y
            %     error('Control matrix not aligned with spatial coordiantes');
            % end
            % this.Bq = Bq;
            this.beta_reg = beta_reg;
            % Quadrature weights for the time integral in the control regularization.
            w = ones(n_t - 1, 1);
            w(1) = 0.5;
            w(2) = 0.5;
            w = (T - this.t_mesh(2)) * w / (n_t - 2);
            this.w_z = repelem(w, n_q);
        end

    end
end
