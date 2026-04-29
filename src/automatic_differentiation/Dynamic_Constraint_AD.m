%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Dynamic_Constraint_AD < Dynamic_Constraint
    % Define constraint equations through an ordinary differential equation
    % dy/dt = f(y(t),z,t),
    %  y(0) = h(z).

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
        path_name
        example_path
    end

    methods (Abstract, Access = public)

        [f] = f_AD(this, y, z, t)

        [h] = h_AD(this, z)

    end

    methods (Access = public)

        function [time_index] = Update_Jacobian(this, y, z, t)
            time_index = find(this.t_mesh - t == 0);
            if ~isempty(time_index)
                if (norm(z - this.z_current(:, time_index)) ~= 0) || (norm(y - this.yt_current(:, time_index)) ~= 0)
                    this.yt_current(:, time_index) = y;
                    this.z_current(:, time_index) = z;
                    [this.Jac_current(:, :, time_index), this.f_current(:, time_index)] = Jac_f_AD_Jac(this, [y; z], t);
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
                    this.Hess_current(:, :, time_index) = Hess_f_AD_Hes(this, [y; z], t, lambda);
                end
            end
        end

        function [f, f_y, f_z] = f(this, y, z, t)
            time_index = this.Update_Jacobian(y, z, t);
            if isempty(time_index)
                [Jac, Fun] = Jac_f_AD_Jac(this, [y; z], t);
                f = Fun;
                f_y = Jac(:, 1:this.n_y);
                f_z = Jac(:, (this.n_y + 1):end);
            else
                f = this.f_current(:, time_index);
                f_y = this.Jac_current(:, 1:this.n_y, time_index);
                f_z = this.Jac_current(:, (this.n_y + 1):end, time_index);
            end
        end

        function [h, h_z] = h(this, z)
            if this.ic_Jac_zero
                h = this.h_AD(z);
                h_z = zeros(length(h), length(z));
            else
                [h_z, h] = Jac_h_AD_Jac(this, z);
            end
        end

        function [Mv] = f_yy_Apply(this, v, y, z, t, lambda)
            time_index = this.Update_Hessian(y, z, t, lambda);
            if isempty(time_index)
                M = Hess_f_AD_Hes(this, [y; z], t, lambda);
                Mv = M(1:this.n_y, 1:this.n_y) * v;
            else
                Mv = this.Hess_current(1:this.n_y, 1:this.n_y, time_index) * v;
            end
        end

        function [Mv] = f_yz_Apply(this, v, y, z, t, lambda)
            time_index = this.Update_Hessian(y, z, t, lambda);
            if isempty(time_index)
                M = Hess_f_AD_Hes(this, [y; z], t, lambda);
                Mv = M(1:this.n_y, (this.n_y + 1):end) * v;
            else
                Mv = this.Hess_current(1:this.n_y, (this.n_y + 1):end, time_index) * v;
            end
        end

        function [Mv] = f_zy_Apply(this, v, y, z, t, lambda)
            time_index = this.Update_Hessian(y, z, t, lambda);
            if isempty(time_index)
                M = Hess_f_AD_Hes(this, [y; z], t, lambda);
                Mv = M((this.n_y + 1):end, 1:this.n_y) * v;
            else
                Mv = this.Hess_current((this.n_y + 1):end, 1:this.n_y, time_index) * v;
            end
        end

        function [Mv] = f_zz_Apply(this, v, y, z, t, lambda)
            time_index = this.Update_Hessian(y, z, t, lambda);
            if isempty(time_index)
                M = Hess_f_AD_Hes(this, [y; z], t, lambda);
                Mv = M((this.n_y + 1):end, (this.n_y + 1):end) * v;
            else
                Mv = this.Hess_current((this.n_y + 1):end, (this.n_y + 1):end, time_index) * v;
            end
        end

        function [Mv] = h_zz_Apply(this, v, z, lambda)
            if this.ic_Hess_zero
                Mv = 0 * v;
            else
                M = Hess_h_AD_Hess(this, z);
                Mv = M * v;
            end
        end

    end

    methods (Access = public)

        function this = Dynamic_Constraint_AD(n_y, n_z, T, n_t)
            arguments
                n_y int32
                n_z int32
                T double
                n_t int32
            end
            this@Dynamic_Constraint(n_y, n_z, T, n_t);

            this.yt_current = inf * ones(n_y, n_t);
            this.z_current = inf * ones(n_z, n_t);
            this.lambdat_current = inf * ones(n_y, n_t);
            this.f_current = inf * ones(n_y, n_t);
            this.Jac_current = inf * ones(n_y, n_y + n_z, n_t);
            this.Hess_current = inf * ones(n_y + n_z, n_y + n_z, n_t);
        end

        function [] = Clear_AD(this)
            evalc("rmpath(this.path_name)");
            evalc("rmpath(this.example_path)");
            evalc("rmdir(this.path_name,'s')");
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

            this.example_path = pwd();
            addpath(this.example_path);
            addpath(this.path_name);
            cd(this.path_name);

            % Test if any functions have changed
            repeat = 1;
            while repeat
                y = randn(this.n_y, 1);
                z = randn(this.n_z, 1);
                t = rand;
                lambda = randn(this.n_y, 1);
                fcurrent = this.f_AD(y, z, t);
                repeat = isnan(norm(fcurrent));
            end
            try
                [~, fold] = Jac_f_AD_Jac(this, [y; z], t);
            catch
                fold = zeros(this.n_y, 1);
            end
            if norm(fold - fcurrent) > 10^-15
                if this.verbose
                    disp('Detected change in constraint');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;

                gyz = adigatorCreateDerivInput([length(y) + length(z), 1], 'yz'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('Jac_f_AD', {this, gyz, t}, options);
                catch
                    if this.verbose
                        disp('RHS jacobian is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('Hess_f_AD', {this, gyz, t, lambda}, options);
                catch
                    if this.verbose
                        disp('RHS Hessian is zero');
                    end
                    this.Hess_zero = true;
                end

            end

            try
                [~, icold] = Jac_h_AD_Jac(this, z);
            catch
                icold = zeros(this.n_y, 1);
            end
            iccurrent = this.h_AD(z);
            if norm(icold - iccurrent) > 10^-15
                if this.verbose
                    disp('Detected change in initial condition');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;
                gz = adigatorCreateDerivInput([length(z), 1], 'z'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('Jac_h_AD', {this, gz}, options);
                catch
                    if this.verbose
                        disp('Initial condition jacobian is zero');
                    end
                    this.ic_Jac_zero = true;
                end

                try
                    genout = adigatorGenHesFile('Hess_h_AD', {this, gz}, options);
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
