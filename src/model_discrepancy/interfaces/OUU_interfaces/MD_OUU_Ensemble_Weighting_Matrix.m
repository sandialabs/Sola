classdef MD_OUU_Ensemble_Weighting_Matrix < handle

    properties
        md_ouu_data_interface
        us_prior_interface
        max_marginal_var_percent
        min_cond_var_percent
        assume_independent

        reg_opt
        reg_opt_min_cond_var_percent

        W_s
        W_s_inv
        R_inv
        C
    end

    methods (Access = public)


        function this = MD_OUU_Ensemble_Weighting_Matrix(md_ouu_data_interface, us_prior_interface, max_marginal_var_percent, min_cond_var_percent, assume_independent)
            arguments
                md_ouu_data_interface MD_OUU_Data_Interface
                us_prior_interface MD_u_Prior_Interface
                max_marginal_var_percent (1, 1) {mustBeNumeric} = 1.0
                min_cond_var_percent (1, 1) {mustBeNumeric} = 0.1
                assume_independent {boolean} = false
            end
            this.md_ouu_data_interface = md_ouu_data_interface;
            this.us_prior_interface = us_prior_interface;
            this.max_marginal_var_percent = max_marginal_var_percent;
            this.min_cond_var_percent = min_cond_var_percent;
            this.assume_independent = assume_independent;

            if this.assume_independent
                ens_size = size(this.md_ouu_data_interface.Reshape_State_to_Mat(this.md_ouu_data_interface.D(:,1)),2);
                this.W_s = eye(ens_size);
                this.W_s_inv = eye(ens_size);
                this.R_inv = eye(ens_size);
                this.C = eye(ens_size);
            else
                this.Compute_Matrices();
            end
        end

        function [] = Compute_Matrices(this)
            ens_size = size(this.md_ouu_data_interface.Reshape_State_to_Mat(this.md_ouu_data_interface.D(:,1)),2);
            N = size(this.md_ouu_data_interface.D,2);
            T = zeros(ens_size,ens_size,N);
            for i = 1:N
                d = this.md_ouu_data_interface.Reshape_State_to_Mat(this.md_ouu_data_interface.D(:,i));
                for s = 1:ens_size
                    for k = 1:ens_size
                        T(s,k,i) = (d(:,s) - d(:,k))' * this.us_prior_interface.Apply_M_u(d(:,s) - d(:,k));
                    end
                    T(s,s,i) = d(:,s)' * this.us_prior_interface.Apply_M_u(d(:,s));
                end
            end
            T = mean(T,3) + 1.e-16;

            A = diag(sqrt(diag(T))) * (1./T) * diag(sqrt(diag(T)));

            [~,i] = min(sum(T,1));
            tmp = T(i,:);
            tmp(i) = Inf;
            [~,j] = min(tmp);

            obj_fun = @(reg_coeff) this.Min_Cond_Var_Obj(A,reg_coeff,i,j);

            x0 = 0.05;
            lb = 0; 
            ub = 1; 
            options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp');
            this.reg_opt = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);

            [this.W_s, this.W_s_inv, this.R_inv, this.C] = this.Assemble_Matrices(A, this.reg_opt);
            this.reg_opt_min_cond_var_percent = (this.W_s_inv(i,i) - this.W_s_inv(i,j)^2/this.W_s_inv(j,j)) / this.W_s_inv(i,i);

            if abs(this.reg_opt_min_cond_var_percent - this.min_cond_var_percent)/abs(this.reg_opt_min_cond_var_percent) > 1.e-2
                disp(['min_cond_var_percent was set to ', num2str(this.min_cond_var_percent), ' but a value of ',num2str(this.reg_opt_min_cond_var_percent),' was achieved'])
            end
        end

        function [Ws, Wsinv, Rinv, C] = Assemble_Matrices(this, A, reg_coeff)
            ens_size = size(A,1);
            C = reg_coeff * (A - eye(ens_size)) + eye(ens_size);
            Ws = diag(diag(C) + 2 * sum(C, 2)) - 2 * C;
            R = chol(Ws);
            Rinv = linsolve(R, eye(ens_size));
            Wsinv = Rinv * Rinv';

            scaling = this.max_marginal_var_percent/max(diag(Wsinv));
            Wsinv = scaling * Wsinv;
            Ws = (1 / scaling) * Ws;
            Rinv = sqrt(scaling) * Rinv;
        end

        function [val] = Min_Cond_Var_Obj(this,A,reg_coeff,i,j)
            [~,Wsinv] = Assemble_Matrices(this, A, reg_coeff);
            val = (Wsinv(i,i) - Wsinv(i,j)^2/Wsinv(j,j) - this.min_cond_var_percent*Wsinv(i,i))^2;
        end

    end

end
