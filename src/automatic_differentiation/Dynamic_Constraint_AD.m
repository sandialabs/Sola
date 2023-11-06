classdef Dynamic_Constraint_AD < Dynamic_Constraint

    properties
        Hess_zero
        ic_Jac_zero
        ic_Hess_zero
    end

    methods (Abstract, Access = public)

        [f] = Time_Instance_RHS_AD(this, y, z, t)

        [h] = Initial_Condition_AD(this, z)

    end

    methods (Access = public)

        function [f, f_y, f_z] = Time_Instance_RHS(this, y, z, t)
            [Jac, Fun] = Jac_Time_Instance_RHS_AD_Jac(this, [y; z], t);

            f = Fun;
            f_y = Jac(:, 1:this.m);
            f_z = Jac(:, (this.m + 1):end);
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
            M = Hess_Time_Instance_RHS_AD_Hes(this, [y; z], t, lambda);
            Mv = M(1:this.m, 1:this.m) * v;
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(this, v, y, z, t, lambda)
            M = Hess_Time_Instance_RHS_AD_Hes(this, [y; z], t, lambda);
            Mv = M(1:this.m, (this.m + 1):end) * v;
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(this, v, y, z, t, lambda)
            M = Hess_Time_Instance_RHS_AD_Hes(this, [y; z], t, lambda);
            Mv = M((this.m + 1):end, 1:this.m) * v;
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(this, v, y, z, t, lambda)
            M = Hess_Time_Instance_RHS_AD_Hes(this, [y; z], t, lambda);
            Mv = M((this.m + 1):end, (this.m + 1):end) * v;
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
    % m: the dimension of the ODE state y(t)
    % n: the dimension of the control z
    % T: the final time
    % N: the number of nodes in the time mesh
    methods (Access = public)

        function this = Dynamic_Constraint_AD(m, n, T, N)
            this@Dynamic_Constraint(m, n, T, N);
        end

        function [] = AD_Initialization(this)
            addpath .;
            if ~exist('AdiGator_Files', 'dir')
                mkdir AdiGator_Files;
            end
            addpath AdiGator_Files;
            cd AdiGator_Files;

            % Test if any functions have changed
            y = randn(this.m, 1);
            z = randn(this.n, 1);
            t = rand;
            lambda = randn(this.m, 1);
            try
                [~, fold] = Jac_Time_Instance_RHS_AD_Jac(this, [y; z], t);
            catch
                fold = zeros(this.m, 1);
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
                icold = zeros(this.m, 1);
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
