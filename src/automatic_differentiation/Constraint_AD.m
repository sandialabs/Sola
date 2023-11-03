classdef Constraint_AD < Constraint

    properties
        n_u
        n_z
        u_current
        c_uu_zero
        c_uz_zero
        c_zz_zero
        verbose
    end

    methods (Abstract, Access = public)

        [c] = constraint_residual(this, u, z)

    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u = this.u_current;
            res = this.constraint_residual(u, z);
            res_norm = norm(res);

            if res_norm > 1.e-12
                J_u = AD_u_constraint_residual_Jac(this, u, z);
                u = u - J_u \ res;
            end

            res = this.constraint_residual(u, z);
            res_norm = norm(res);
            if res_norm > 1.e-12
                % Execute nonlinear solve to determine the state
                options = optimoptions('fsolve', 'Display', 'none', 'OptimalityTolerance', 1.e-14, 'SpecifyObjectiveGradient', true, 'CheckGradients', false);
                u = fsolve(@(u)this.Nonlinear_Solver_Evaluation(u, z), u, options);
            end

            this.u_current = u;
        end

        function [res, J_u] = Nonlinear_Solver_Evaluation(this, u, z)
            [J_u, res] = AD_u_constraint_residual_Jac(this, u, z);
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            J_u = AD_u_constraint_residual_Jac(this, u, z);
            Mv = J_u' \ v;
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            J_z = AD_z_constraint_residual_Jac(this, u, z);
            Mv = J_z' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            J_u = AD_u_constraint_residual_Jac(this, u, z);
            Mv = J_u \ v;
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            J_z = AD_z_constraint_residual_Jac(this, u, z);
            Mv = J_z * v;
        end

        function [c] = c(this, u, z)
            c = this.constraint_residual(u, z);
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            if ~this.c_uu_zero
                M = AD_uu_constraint_residual_Hes(this, u, z, lambda);
                Mv = M * v;
            else
                Mv = 0 * v;
            end
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            if ~this.c_uz_zero
                M = AD_uz_constraint_residual_Hes(this, [u; z], lambda);
                Mv = M(1:this.n_u, (this.n_u + 1):end) * v;
            else
                Mv = 0 * u;
            end
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            if ~this.c_uz_zero
                M = AD_uz_constraint_residual_Hes(this, [u; z], lambda);
                Mv = M((this.n_u + 1):end, 1:this.n_u) * v;
            else
                Mv = 0 * z;
            end
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            if ~this.c_zz_zero
                M = AD_zz_constraint_residual_Hes(this, u, z, lambda);
                Mv = M * v;
            else
                Mv = 0 * v;
            end
        end

    end

    methods (Access = public)

        function this = Constraint_AD(n_u, n_z)
            this.n_u = n_u;
            this.n_z = n_z;
            this.u_current = zeros(this.n_u, 1);
            this.c_uu_zero = false;
            this.c_uz_zero = false;
            this.c_zz_zero = false;
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
                [~, cold] = AD_u_constraint_residual_Jac(this, u, z);
            catch
                cold = zeros(this.n_u, 1);
            end
            ccurrent = this.constraint_residual(u, z);
            if norm(cold - ccurrent) > 10^-15
                if this.verbose
                    disp('Detected change in constraint');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;
                gu = adigatorCreateDerivInput([length(u), 1], 'u'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('AD_u_constraint_residual', {this, gu, z}, options);
                catch
                    if this.verbose
                        disp('constraint u jacobian is zero');
                    end
                end

                gz = adigatorCreateDerivInput([length(z), 1], 'z'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('AD_z_constraint_residual', {this, u, gz}, options);
                catch
                    if this.verbose
                        disp('constraint z jacobian is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('AD_uu_constraint_residual', {this, gu, z, lambda}, options);
                catch
                    if this.verbose
                        disp('constraint uu hessian is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('AD_zz_constraint_residual', {this, u, gz, lambda}, options);
                catch
                    if this.verbose
                        disp('constraint zz hessian is zero');
                    end
                end

                guz = adigatorCreateDerivInput([length(u) + length(z), 1], 'uz'); % Create Deriv Input
                try
                    genout = adigatorGenHesFile('AD_uz_constraint_residual', {this, guz, lambda}, options);
                catch
                    if this.verbose
                        disp('constraint uz hessian is zero');
                    end
                end

            end

            try
                AD_uu_constraint_residual_Hes(this, u, z, lambda);
            catch
                if this.verbose
                    disp('constraint uu hessian is zero');
                end
                this.c_uu_zero = true;
            end

            try
                AD_zz_constraint_residual_Hes(this, u, z, lambda);
            catch
                if this.verbose
                    disp('constraint zz hessian is zero');
                end
                this.c_zz_zero = true;
            end

            try
                AD_uz_constraint_residual_Hes(this, [u; z], lambda);
            catch
                if this.verbose
                    disp('constraint uz hessian is zero');
                end
                this.c_uz_zero = true;
            end

            cd ..;
        end

    end

end
