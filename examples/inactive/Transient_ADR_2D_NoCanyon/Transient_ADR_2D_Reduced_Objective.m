%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Transient_ADR_2D_Reduced_Objective < Dynamic_Objective
    % Objective function for a two-species time-dependent
    % advection-diffusion-reaction problem in two spatial dimensions
    % with separate reduced-order approximations for each state variable.
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
    %
    % The states are approximated linearly with
    %
    % .. math::
    %    \y_1(t) \approx \V_1\yhat_1(t),
    %    \qquad
    %    \y_2(t) \approx \V_2\yhat_2(t),
    %
    % where :math:`\y_i(t) \in \R^{r_i}` and :math:`\V_i \in \R^{n_x \times r_i}`.

    properties (SetAccess = protected)
        fullobj         % full-space objective, a :class:`Transient_ADR_2D_Objective` object.
        r_2             % Dimension :math:`r_2` of the reduced state :math:`\yhat_2`.
        r_1             % Dimension :math:`r_1` of the reduced state :math:`\yhat_1`.
        V1pMpV1         % reduced weighted mass matrix :math:`\V_1\trp(\p\trp\ast\M\ast\p)\V_1`.
    end

    methods (Access = public)

        function [val, grad_y] = g(this, y, ~)
            % :math:`g(yhat) = || V1 yhat_1(t) .* p ||_M^2 / 2`
            % :math:`g_y(yhat) = V1pMpV1 yhat_1(t)

            y1 = reshape(y(1:this.r_1), this.r_1, 1);
            grad_y1 = this.V1pMpV1 * y1;
            grad_y2 = zeros(this.r_2, 1);
            val = .5 * y1' * grad_y1;
            grad_y = [grad_y1; grad_y2];
        end

        function [val, grad_z] = R(this, z)
            [val, grad_z] = this.fullobj.R(z);
        end

        function [y_out] = g_yy_Apply(this, y_in, ~, ~)
            y_in1 = y_in(1:this.r_1, :);
            y_out = zeros(this.n_y, size(y_in, 2));
            y_out(1:this.r_1, :) = this.V1pMpV1 * y_in1;
        end

        function [z_out] = R_zz_Apply(this, z_in, z)
            z_out = this.fullobj.R_zz_Apply(z_in, z);
        end

    end

    methods

        function this = Transient_ADR_2D_Reduced_Objective(fullobj, V1, V2)
            % Constructor.
            %
            % Parameters
            % ----------
            % objective
            %   Initialized objective for the full-order states, of type
            %   :class:`Transient_ADR_2D_Objective`.
            % V1
            %   Basis matrix :math:`\V_1\in\R^{n_x \times r_1}`.
            % V2
            %   Basis matrix :math:`\V_2\in\R^{n_x \times r_2}`.

            r1 = size(V1, 2);
            r2 = size(V2, 2);
            r = r1 + r2;

            this@Dynamic_Objective(r, fullobj.n_z, fullobj.T, fullobj.n_t);

            this.fullobj = fullobj;
            this.r_1 = r1;
            this.r_2 = r2;
            this.V1pMpV1 = V1' * fullobj.pMp * V1;
        end

    end
end
