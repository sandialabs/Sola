classdef Dynamic_Objective_AD < Dynamic_Objective
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
    %
    % The gradients and Hessian actions of :math:`g` are computed via
    % automatic differentiation from the :meth:`Time_Instance_Objective_AD()`
    % and :meth:`Regularization_Objective_AD()` methods.

    properties
        verbose
        path_name
    end

    methods (Abstract, Access = public)

        [val] = Time_Instance_Objective_AD(this, y, t)
        % Evaluate the integrand :math:`g(\y,t)`.
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

        [val] = Regularization_Objective_AD(this, z)
        % Evaluate the regularization term :math:`R(\z)`.
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

    end

    methods (Access = public)

        function [val, grad_y] = Time_Instance_Objective(this, y, t)
            [grad_y, val] = grad_Time_Instance_Objective_AD_Jac(this, y, t);
            grad_y = grad_y';
        end

        function [val, grad_z] = Regularization_Objective(this, z)
            [grad_z, val] = grad_Regularization_Objective_AD_Jac(this, z);
            grad_z = grad_z';
        end

        function [Mv] = Time_Instance_Objective_yy_Apply(this, v, y, t)
            M = Hess_Time_Instance_Objective_AD_Hes(this, y, t);
            Mv = M * v;
        end

        function [Mv] = Regularization_Objective_zz_Apply(this, v, z)
            M = Hess_Regularization_Objective_AD_Hes(this, z);
            Mv = M * v;
        end

    end

    methods (Access = public)

        function this = Dynamic_Objective_AD(n_y, n_z, T, n_t)
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
            this@Dynamic_Objective(n_y, n_z, T, n_t);
            this.verbose = true;
        end

        function [] = AD_Initialization(this, folder_name)

            if nargin > 1
                this.path_name = [pwd(), '/', folder_name];
            else
                this.path_name = [pwd(), '/AdiGator_Files'];
            end

            if ~exist(this.path_name, 'dir')
                mkdir(this.path_name);
            end

            addpath('.');
            addpath(this.path_name);
            cd(this.path_name);

            % Test if any functions have changed
            y = randn(this.n_y, 1);
            z = randn(this.n_z, 1);
            t = rand;
            try
                [~, valold] = grad_Time_Instance_Objective_AD_Jac(this, y, t);
            catch
                valold = 0;
            end
            valcurrent = this.Time_Instance_Objective_AD(y, t);
            if norm(valold - valcurrent) > 10^-15
                if this.verbose
                    disp('Detected change in time instance objective');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;
                gy = adigatorCreateDerivInput([length(y), 1], 'y'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('grad_Time_Instance_Objective_AD', {this, gy, t}, options);
                catch
                    if this.verbose
                        disp('objective gradient is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('Hess_Time_Instance_Objective_AD', {this, gy, t}, options);
                catch
                    if this.verbose
                        disp('objective Hessian is zero');
                    end
                end

            end

            try
                [~, valold] = grad_Regularization_Objective_AD_Jac(this, z);
            catch
                valold = 0;
            end
            valcurrent = this.Regularization_Objective_AD(z);
            if norm(valold - valcurrent) > 10^-15
                if this.verbose
                    disp('Detected change in regularization objective');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;
                gz = adigatorCreateDerivInput([length(z), 1], 'z'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('grad_Regularization_Objective_AD', {this, gz}, options);
                catch
                    if this.verbose
                        disp('regularization gradient is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('Hess_Regularization_Objective_AD', {this, gz}, options);
                catch
                    if this.verbose
                        disp('regularization Hessian is zero');
                    end
                end

            end
            cd ..;
        end

    end

end
