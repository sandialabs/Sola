classdef Objective_AD < Objective

    properties
        n_u
        n_z
        verbose
    end

    methods (Abstract, Access = public)

        [val] = J_val(this, u, z)

    end

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            [grad, val] = AD_grad_J_val_Jac(this, [u; z]);
            grad_u = grad(1:this.n_u)';
            grad_z = grad((this.n_u + 1):end)';
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            H = AD_Hess_J_val_Hes(this, [u; z]);
            Mv = H(1:this.n_u, 1:this.n_u) * v;
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            H = AD_Hess_J_val_Hes(this, [u; z]);
            Mv = H(1:this.n_u, (this.n_u + 1):end) * v;
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            H = AD_Hess_J_val_Hes(this, [u; z]);
            Mv = H((this.n_u + 1):end, 1:this.n_u) * v;
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            H = AD_Hess_J_val_Hes(this, [u; z]);
            Mv = H((this.n_u + 1):end, (this.n_u + 1):end) * v;
        end

    end

    methods (Access = public)

        function this = Objective_AD(n_u, n_z)
            this.n_u = n_u;
            this.n_z = n_z;
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
            try
                [~, valold] = AD_grad_J_val_Jac(this, [u; z]);
            catch
                valold = 0;
            end
            valcurrent = this.J_val(u, z);
            if norm(valold - valcurrent) > 10^-15
                if this.verbose
                    disp('Detected change in objective');
                end

                options = adigatorOptions();
                options.overwrite = 1;
                options.echo = 0;
                guz = adigatorCreateDerivInput([length(u) + length(z), 1], 'uz'); % Create Deriv Input
                try
                    genout = adigatorGenJacFile('AD_grad_J_val', {this, guz}, options);
                catch
                    if this.verbose
                        disp('objective gradient is zero');
                    end
                end

                try
                    genout = adigatorGenHesFile('AD_Hess_J_val', {this, guz}, options);
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
