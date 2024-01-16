classdef Reduced_Dynamic_Objective < Dynamic_Objective
    % For an objective function as defined by :class:`Dynamic_Objective`,
    %
    % .. math:: J(\u,\z) = \sum_{j=1}^{N} w_{j} g(\y_{j}, t_{j}) + R(\z),
    %
    % this class defines a modified objective that results from approximating
    % the ODE state as :math:`\y(t) \approx \V\hat{\y}(t)`:
    %
    % .. math::
    %  \tilde{J}(\hat{\u},\z) = \sum_{j=1}^{N} w_{j} \hat{g}(\hat{\y}_{j}, t_{j}) + R(\z),
    %  \qquad
    %  \hat{g}(\hat{\y}, t) = g(\V\hat{\y}, t),
    %
    % where
    %
    % * :math:`\hat{\u} = (\hat{\y}_{1}\trp \cdots \hat{\y}_{n_t}\trp)\trp\in\R^{n_u'}`
    %   with :math:`\hat{\y}_{j} \in \R^{n_y'}` and :math:`n_u' = n_y' n_t`,
    % * :math:`\V\in\R^{n_y \times n_y'}`,
    % * :math:`0 = t_{1} < t_{2} < \cdots < t_{n_t} = T`
    %   with equal spacing :math:`\delta t = t_{j+1} - t_{j}`,
    % * :math:`\z \in \R^{n_z}`,
    % * :math:`\hat{g}:\R^{n_y'} \times [0, T] \to \R`,
    % * :math:`R:\R^{n_z} \to \R`, and
    % * :math:`w_{j} = \delta t` for :math:`j = 2, \ldots, N - 1`
    %   and :math:`w_{1} = w_{n_t} = \frac{1}{2}\delta t`.

    properties (Access = protected)
        obj
        V
    end

    methods (Access = public)

        %% Constructor.

        function this = Reduced_Dynamic_Objective(objective, V)
            % Parameters
            % ----------
            % objective
            %   Initialized dynamic objective, subclassed from
            %   :class:`Dynamic_Objective`.
            % V
            %   Basis matrix :math:`\V\in\R^{n_y \times n_y'}`.
            this@Dynamic_Objective(objective.n_y, objective.n_z, objective.T, objective.n_t);
            this.obj = objective;
            this.V = V;
        end

        %% Implement abstract methods from parent class.

        function [val, grad_yhat] = Time_Instance_Objective(this, yhat, t)
            % Evaluate the integrand :math:`\hat{g}(\hat{\y},t) = g(\V\hat{\y}, t)`
            % and its gradient
            % :math:`\grad{\hat{y}}\hat{g}(\hat{\y},t) = \V\trp\grad{y}g(\V\hat{\y},t)`.
            %
            % Parameters
            % ----------
            % yhat
            %   Reduced differential equation state :math:`\hat{\y}\in\R^{n_y'}`.
            % t
            %   Time :math:`t`.
            %
            % Returns
            % -------
            % val : double
            %   Function value :math:`g(\y,t)\in\R`.
            % grad_y : vector
            %   Function gradient :math:`\grad{y}g(\y,t)\in\R^{n_y}`.
            [val, grad_y] = this.obj.Time_Instance_Objective(this.V * yhat, t);
            grad_yhat = this.V' * grad_y;
        end

        function [val, grad_z] = Regularization_Objective(this, z)
            [val, grad_z] = this.obj.Regularization_Objective(z);
        end

        function [Mv] = Time_Instance_Objective_yy_Apply(this, vhat, yhat, t)
            % Compute the Hessian-vector product
            % :math:`\grad{\hat{y},\hat{y}}\hat{g}(\hat{\y}, t)\hat{\v}
            % = \V\trp \grad{y,y}g(\V\hat{\y}, t)\V\hat{\v}`.
            %
            % Parameters
            % ----------
            % vhat
            %   Search direction :math:`\hat{\v}\in\R^{n_y'}`.
            % yhat
            %   Reduced differential equation state :math:`\hat{\y}\in\R^{n_y'}`.
            % t
            %   Time :math:`t`.
            %
            % Returns
            % -------
            % Mv : vector
            %   Hessian-vector product
            %   :math:`\grad{\hat{y},\hat{y}}\hat{g}(\hat{\y}, t)\hat{\v}\in\R^{n_y'}`.
            v = this.V * vhat;
            y = this.V * yhat;
            Mv = this.V' * this.obj.Time_Instance_Objective_yy_Apply(v, y, t);
        end

        function [Mv] = Regularization_Objective_zz_Apply(this, v, z)
            Mv = this.obj.Regularization_Objective_zz_Apply(v, z);
        end

    end
end
