classdef Dynamic_Objective_AD < Dynamic_Objective

    properties
        verbose
    end

    methods (Abstract, Access = public)

        [val] = Time_Instance_Objective_AD(this, y, t)

        [val] = Regularization_Objective_AD(this, z)

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

        % Input:
        % n_y: the dimension of the ODE state y(t)
        % n_z: the dimension of the control z
        % T: the final time
        % n_t: the number of nodes in the time mesh
        function this = Dynamic_Objective_AD(n_y, n_z, T, n_t)
            this@Dynamic_Objective(n_y, n_z, T, n_t);
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
