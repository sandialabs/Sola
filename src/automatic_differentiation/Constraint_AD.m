classdef Constraint_AD < Constraint
    % Define a constraint function :math:`\c(\u, \z)\to\R^{n_u}` where
    % :math:`\u \in \R^{n_u}` is the state and
    % :math:`\z \in \R^{n_z}` is the control.
    %
    % The Jacobian and Hessian actions of $\c$ are computed via
    % automatic differentiation from the :meth:`c_AD()` method.

    properties
        n_u             % Dimension :math:`n_u` of the state :math:`\u`.
        n_z             % Dimension :math:`n_z` of the control :math:`\z`.
        u_current
        z_current
        lambda_current
        Jac_current
        Hess_current
        Hess_zero
        verbose         % If ``true``, print automatic differentiation info.
    end

    methods (Abstract, Access = public)

        [c] = c_AD(this, u, z)
        % Constraint :math:`\c(\u,\z)`.
        % This method is used for automatic differentiation.
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
        % Mv : vector
        %   Constraint :math:`\c(\u,\z)\in\R^{n_u}`.

    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            % Given :math:`\z`, solve the constraint equation
            % :math:`\c(\u,\z)=\0` for :math:`\u`, i.e., compute
            % :math:`\u = \S(\z)`.
            %
            % Parameters
            % ----------
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            %
            % Returns
            % -------
            % u : vector
            %   State :math:`\u = \S(\z) \in \R^{n_u}`.

            u = this.u_current;
            res = this.c_AD(u, z);

            if norm(res) > 1.e-12
                J = Jac_c_AD_Jac(this, [u; z]);
                J_u = J(:, 1:this.n_u);
                u = u - J_u \ res;

                res = this.c_AD(u, z);
                res_norm = norm(res);
                if res_norm > 1.e-12
                    % Execute nonlinear solve to determine the state
                    options = optimoptions('fsolve', 'Display', 'none', 'OptimalityTolerance', 1.e-14, 'SpecifyObjectiveGradient', true, 'CheckGradients', false);
                    u = fsolve(@(u)this.Nonlinear_Solver_Evaluation(u, z), u, options);
                end
            end

        end

        function [] = Update_Jacobian(this, u, z)
            if (norm(z - this.z_current) ~= 0) || (norm(u - this.u_current) ~= 0)
                this.u_current = u;
                this.z_current = z;
                this.Jac_current = Jac_c_AD_Jac(this, [u; z]);
            end
        end

        function [] = Update_Hessian(this, u, z, lambda)
            if (norm(z - this.z_current) ~= 0) || (norm(u - this.u_current) ~= 0) || (norm(lambda - this.lambda_current) ~= 0)
                this.u_current = u;
                this.z_current = z;
                this.lambda_current = lambda;
                this.Hess_current = Hess_c_AD_Hes(this, [u; z], lambda);
            end
        end

        function [res, J_u] = Nonlinear_Solver_Evaluation(this, u, z)
            [J, res] = Jac_c_AD_Jac(this, [u; z]);
            J_u = J(:, 1:this.n_u);
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            this.Update_Jacobian(u, z);
            Mv = this.Jac_current(:, 1:this.n_u)' \ v;
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            this.Update_Jacobian(u, z);
            Mv = this.Jac_current(:, (this.n_u + 1):end)' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            this.Update_Jacobian(u, z);
            Mv = this.Jac_current(:, 1:this.n_u) \ v;
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            this.Update_Jacobian(u, z);
            Mv = this.Jac_current(:, (this.n_u + 1):end) * v;
        end

        function [c] = c(this, u, z)
            c = this.c_AD(u, z);
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            if ~this.Hess_zero
                this.Update_Hessian(u, z, lambda);
                Mv = this.Hess_current(1:this.n_u, 1:this.n_u) * v;
            else
                Mv = 0 * v;
            end
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            if ~this.Hess_zero
                this.Update_Hessian(u, z, lambda);
                Mv = this.Hess_current(1:this.n_u, (this.n_u + 1):end) * v;
            else
                Mv = 0 * u;
            end
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            if ~this.Hess_zero
                this.Update_Hessian(u, z, lambda);
                Mv = this.Hess_current((this.n_u + 1):end, 1:this.n_u) * v;
            else
                Mv = 0 * z;
            end
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            if ~this.Hess_zero
                this.Update_Hessian(u, z, lambda);
                Mv = this.Hess_current((this.n_u + 1):end, (this.n_u + 1):end) * v;
            else
                Mv = 0 * v;
            end
        end

    end

    methods (Access = public)

        function this = Constraint_AD(n_u, n_z)
            % Parameters
            % ----------
            % n_u
            %   Dimension :math:`n_u` of the state.
            % n_z
            %   Dimension :math:`n_z` of the control.
            this.n_u = n_u;
            this.n_z = n_z;
            this.u_current = inf * ones(this.n_u, 1);
            this.z_current = inf * ones(this.n_z, 1);
            this.lambda_current = inf * ones(this.n_u, 1);
            this.Jac_current = inf * ones(this.n_u, this.n_u + this.n_z);
            this.Hess_current = inf * ones(this.n_u + this.n_z, this.n_u + this.n_z);
            this.Hess_zero = false;
            this.verbose = true;
        end

        function [] = AD_Initialization(this)
            addpath .;
            if ~exist('AdiGator_Files', 'dir')
                mkdir AdiGator_Files;
            end
            addpath AdiGator_Files;
            cd AdiGator_Files;

            % Test if any functions have changed
            u = randn(this.n_u, 1);
            z = randn(this.n_z, 1);
            lambda = randn(this.n_u, 1);
            try
                [~, cold] = Jac_c_AD_Jac(this, [u; z]);
            catch
                cold = zeros(this.n_u, 1);
            end
            ccurrent = this.c_AD(u, z);
            if norm(cold - ccurrent) > 10^-15
                if this.verbose
                    disp('Detected change in constraint');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;
                guz = adigatorCreateDerivInput([length(u) + length(z), 1], 'uz'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('Jac_c_AD', {this, guz}, options);
                catch
                    if this.verbose
                        disp('constraint jacobian is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('Hess_c_AD', {this, guz, lambda}, options);
                catch
                    if this.verbose
                        disp('constraint hessian is zero');
                    end
                    this.Hess_zero = true;
                end

            end

            cd ..;
        end

    end

end
