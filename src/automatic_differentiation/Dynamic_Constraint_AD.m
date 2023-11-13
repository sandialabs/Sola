classdef Dynamic_Constraint_AD < Dynamic_Constraint

    properties
        yt_current
        z_current
        lambdat_current
        f_current
        Jac_current
        Hess_current
        Hess_zero
        ic_Jac_zero
        ic_Hess_zero
    end

    methods (Abstract, Access = public)

        [f] = Time_Instance_RHS_AD(this, y, z, t)

        [h] = Initial_Condition_AD(this, z)

    end

    methods (Access = public)

        function [time_index] = Update_Jacobian(this, y, z, t)
            time_index = find(this.t_mesh - t == 0);
            if ~isempty(time_index)
                if (norm(z - this.z_current(:, time_index)) ~= 0) || (norm(y - this.yt_current(:, time_index)) ~= 0)
                    this.yt_current(:, time_index) = y;
                    this.z_current(:, time_index) = z;
                    [this.Jac_current(:, :, time_index), this.f_current(:, time_index)] = Jac_Time_Instance_RHS_AD_Jac(this, [y; z], t);
                end
            end
        end

        function [time_index] = Update_Hessian(this, y, z, t, lambda)
            time_index = find(this.t_mesh - t == 0);
            if ~isempty(time_index)
                if (norm(z - this.z_current(:, time_index)) ~= 0) || (norm(y - this.yt_current(:, time_index)) ~= 0) || (norm(lambda - this.lambdat_current(:, time_index)) ~= 0)
                    this.yt_current(:, time_index) = y;
                    this.z_current(:, time_index) = z;
                    this.lambdat_current(:, time_index) = lambda;
                    this.Hess_current(:, :, time_index) = Hess_Time_Instance_RHS_AD_Hes(this, [y; z], t, lambda);
                end
            end
        end

        function [f, f_y, f_z] = Time_Instance_RHS(this, y, z, t)
            time_index = this.Update_Jacobian(y, z, t);
            if isempty(time_index)
                [Jac, Fun] = Jac_Time_Instance_RHS_AD_Jac(this, [y; z], t);
                f = Fun;
                f_y = Jac(:, 1:this.n_y);
                f_z = Jac(:, (this.n_y + 1):end);
            else
                f = this.f_current(:, time_index);
                f_y = this.Jac_current(:, 1:this.n_y, time_index);
                f_z = this.Jac_current(:, (this.n_y + 1):end, time_index);
            end
        end

        function [h, h_z] = Initial_Condition(this, z)
            if this.ic_Jac_zero
                h = this.Initial_Condition_AD(z);
                h_z = zeros(length(h), length(z));
            else
                [h_z, h] = Jac_Initial_Condition_AD_Jac(this, z);
            end
        end

        function [Mv] = Time_Instance_RHS_yy_Apply(this, v, y, z, t, lambda)
            time_index = this.Update_Hessian(y, z, t, lambda);
            if isempty(time_index)
                M = Hess_Time_Instance_RHS_AD_Hes(this, [y; z], t, lambda);
                Mv = M(1:this.n_y, 1:this.n_y) * v;
            else
                Mv = this.Hess_current(1:this.n_y, 1:this.n_y, time_index) * v;
            end
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(this, v, y, z, t, lambda)
            time_index = this.Update_Hessian(y, z, t, lambda);
            if isempty(time_index)
                M = Hess_Time_Instance_RHS_AD_Hes(this, [y; z], t, lambda);
                Mv = M(1:this.n_y, (this.n_y + 1):end) * v;
            else
                Mv = this.Hess_current(1:this.n_y, (this.n_y + 1):end, time_index) * v;
            end
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(this, v, y, z, t, lambda)
            time_index = this.Update_Hessian(y, z, t, lambda);
            if isempty(time_index)
                M = Hess_Time_Instance_RHS_AD_Hes(this, [y; z], t, lambda);
                Mv = M((this.n_y + 1):end, 1:this.n_y) * v;
            else
                Mv = this.Hess_current((this.n_y + 1):end, 1:this.n_y, time_index) * v;
            end
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(this, v, y, z, t, lambda)
            time_index = this.Update_Hessian(y, z, t, lambda);
            if isempty(time_index)
                M = Hess_Time_Instance_RHS_AD_Hes(this, [y; z], t, lambda);
                Mv = M((this.n_y + 1):end, (this.n_y + 1):end) * v;
            else
                Mv = this.Hess_current((this.n_y + 1):end, (this.n_y + 1):end, time_index) * v;
            end
        end

        function [Mv] = Initial_Condition_zz_Apply(this, v, z, lambda)
            if this.ic_Hess_zero
                Mv = 0 * v;
            else
                M = Hess_Initial_Condition_AD_Hess(this, z);
                Mv = M * v;
            end
        end

    end

    % Input:
    % n_y: the dimension of the ODE state y(t)
    % n_z: the dimension of the control z
    % T: the final time
    % n_t: the number of nodes in the time mesh
    methods (Access = public)

        function this = Dynamic_Constraint_AD(n_y, n_z, T, n_t)
            this@Dynamic_Constraint(n_y, n_z, T, n_t);

            this.yt_current = inf * ones(n_y, n_t);
            this.z_current = inf * ones(n_z, n_t);
            this.lambdat_current = inf * ones(n_y, n_t);
            this.f_current = inf * ones(n_y, n_t);
            this.Jac_current = inf * ones(n_y, n_y + n_z, n_t);
            this.Hess_current = inf * ones(n_y + n_z, n_y + n_z, n_t);
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
            lambda = randn(this.n_y, 1);
            try
                [~, fold] = Jac_Time_Instance_RHS_AD_Jac(this, [y; z], t);
            catch
                fold = zeros(this.n_y, 1);
            end
            fcurrent = this.Time_Instance_RHS_AD(y, z, t);
            if norm(fold - fcurrent) > 10^-15
                if this.verbose
                    disp('Detected change in constraint');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;

                gyz = adigatorCreateDerivInput([length(y) + length(z), 1], 'yz'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('Jac_Time_Instance_RHS_AD', {this, gyz, t}, options);
                catch
                    if this.verbose
                        disp('RHS jacobian is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('Hess_Time_Instance_RHS_AD', {this, gyz, t, lambda}, options);
                catch
                    if this.verbose
                        disp('RHS Hessian is zero');
                    end
                    this.Hess_zero = true;
                end

            end

            try
                [~, icold] = Jac_Initial_Condition_AD_Jac(this, z);
            catch
                icold = zeros(this.n_y, 1);
            end
            iccurrent = this.Initial_Condition_AD(z);
            if norm(icold - iccurrent) > 10^-15
                if this.verbose
                    disp('Detected change in initial condition');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;
                gz = adigatorCreateDerivInput([length(z), 1], 'z'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('Jac_Initial_Condition_AD', {this, gz}, options);
                catch
                    if this.verbose
                        disp('Initial condition jacobian is zero');
                    end
                    this.ic_Jac_zero = true;
                end

                try
                    genout = adigatorGenHesFile('Hess_Initial_Condition_AD', {this, gz}, options);
                catch
                    if this.verbose
                        disp('Initial condition hessian is zero');
                    end
                    this.ic_Hess_zero = true;
                end
            end

            cd ..;
        end

    end

end
