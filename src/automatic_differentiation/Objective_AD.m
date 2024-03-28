classdef Objective_AD < Objective
    % Define a scalar-valued objective function
    %
    % .. math:: J(\u, \z) \to \R
    %
    % where
    % :math:`\u \in \R^{n_u}` is the state and
    % :math:`\z \in \R^{n_z}` is the control.
    %
    % The gradients and Hessian actions of :math:`J` are computed via
    % automatic differentiation from the :meth:`J_AD()` method.

    properties
        n_u             % Dimension :math:`n_u` of the state :math:`\u`.
        n_z             % Dimension :math:`n_z` of the control :math:`\z`.
        verbose         % If ``true``, print automatic differentiation info.
        path_name
    end

    methods (Abstract, Access = public)

        [val] = J_AD(this, u, z)
        % Evaluate the objective function :math:`J(\u,\z)`.
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
        % val : double
        %   Objective value :math:`J(\u,\z)`.

    end

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            [grad, val] = grad_J_AD_Jac(this, [u; z]);
            grad_u = grad(1:this.n_u)';
            grad_z = grad((this.n_u + 1):end)';
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            H = Hess_J_AD_Hes(this, [u; z]);
            Mv = H(1:this.n_u, 1:this.n_u) * v;
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            H = Hess_J_AD_Hes(this, [u; z]);
            Mv = H(1:this.n_u, (this.n_u + 1):end) * v;
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            H = Hess_J_AD_Hes(this, [u; z]);
            Mv = H((this.n_u + 1):end, 1:this.n_u) * v;
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            H = Hess_J_AD_Hes(this, [u; z]);
            Mv = H((this.n_u + 1):end, (this.n_u + 1):end) * v;
        end

    end

    methods (Access = public)

        function this = Objective_AD(n_u, n_z)
            this.n_u = n_u;
            this.n_z = n_z;
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
            u = randn(this.n_u, 1);
            z = randn(this.n_z, 1);
            try
                [~, valold] = grad_J_AD_Jac(this, [u; z]);
            catch
                valold = 0;
            end
            valcurrent = this.J_AD(u, z);
            if norm(valold - valcurrent) > 10^-15
                if this.verbose
                    disp('Detected change in objective');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;
                guz = adigatorCreateDerivInput([length(u) + length(z), 1], 'uz'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('grad_J_AD', {this, guz}, options);
                catch
                    if this.verbose
                        disp('objective gradient is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('Hess_J_AD', {this, guz}, options);
                catch
                    if this.verbose
                        disp('objective Hessian is zero');
                    end
                end

            end
            cd ..;
        end

    end
end
