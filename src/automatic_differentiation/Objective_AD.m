%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Objective_AD < Objective
    % Define a scalar-valued objective function J(u,z)

    properties
        n_u             % Dimension of the state u.
        n_z             % Dimension of the control z.
        verbose         % If ``true``, print automatic differentiation info.
        path_name
    end

    methods (Abstract, Access = public)

        [val] = J_AD(this, u, z)

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
