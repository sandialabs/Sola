%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Dynamic_Objective_AD < Dynamic_Objective
    % Define an objective function J(u,z) = \int_{0}^{T} g(y(t), t) dt + R(z).

    properties
        verbose
        path_name
    end

    methods (Abstract, Access = public)

        [val] = g_AD(this, y, t)

        [val] = R_AD(this, z)

    end

    methods (Access = public)

        function [val, grad_y] = g(this, y, t)
            [grad_y, val] = grad_g_AD_Jac(this, y, t);
            grad_y = grad_y';
        end

        function [val, grad_z] = R(this, z)
            [grad_z, val] = grad_R_AD_Jac(this, z);
            grad_z = grad_z';
        end

        function [Mv] = g_yy_Apply(this, v, y, t)
            M = Hess_g_AD_Hes(this, y, t);
            Mv = M * v;
        end

        function [Mv] = R_zz_Apply(this, v, z)
            M = Hess_R_AD_Hes(this, z);
            Mv = M * v;
        end

    end

    methods (Access = public)

        function this = Dynamic_Objective_AD(n_y, n_z, T, n_t)
            arguments
                n_y int32
                n_z int32
                T double
                n_t int32
            end
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
                [~, valold] = grad_g_AD_Jac(this, y, t);
            catch
                valold = 0;
            end
            valcurrent = this.g_AD(y, t);
            if norm(valold - valcurrent) > 10^-15
                if this.verbose
                    disp('Detected change in time instance objective');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;
                gy = adigatorCreateDerivInput([length(y), 1], 'y'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('grad_g_AD', {this, gy, t}, options);
                catch
                    if this.verbose
                        disp('objective gradient is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('Hess_g_AD', {this, gy, t}, options);
                catch
                    if this.verbose
                        disp('objective Hessian is zero');
                    end
                end

            end

            try
                [~, valold] = grad_R_AD_Jac(this, z);
            catch
                valold = 0;
            end
            valcurrent = this.R_AD(z);
            if norm(valold - valcurrent) > 10^-15
                if this.verbose
                    disp('Detected change in regularization objective');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;
                gz = adigatorCreateDerivInput([length(z), 1], 'z'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('grad_R_AD', {this, gz}, options);
                catch
                    if this.verbose
                        disp('regularization gradient is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('Hess_R_AD', {this, gz}, options);
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
