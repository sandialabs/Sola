classdef MD_Numeric_Laplacian_u_Prior_Interface < MD_Elliptic_u_Prior_Interface

    properties
        M
        S
        u_hyperparam_interface
        determine_u_hyperparams
        beta_u
        E_u

        is_sparse
        R
        P
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            if this.is_sparse
                tmp1 = this.R' \ (this.P' * u_in);
                tmp2 = this.R \ tmp1;
                u_out = this.P * tmp2;
            else
                tmp = this.R' \ u_in;
                u_out = this.R \ tmp;
            end
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            if this.is_sparse
                tmp1 = this.R' \ (this.P' * u_in);
                tmp2 = this.R \ tmp1;
                u_out = this.P * tmp2;
            else
                tmp = this.R' \ u_in;
                u_out = this.R \ tmp;
            end
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [] = Set_beta_u(this, beta_u_new)
            this.beta_u = beta_u_new;
            this.E_u = this.beta_u * this.S + this.M;
            this.is_sparse = issparse(this.E_u);
            if this.is_sparse
                [this.R, flag, this.P] = chol(this.E_u);
                if flag ~= 0
                    disp('Error in Cholesky factorization of E_u');
                end
            else
                this.R = chol(this.E_u);
            end
        end

        function this = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface)
            this@MD_Elliptic_u_Prior_Interface(u_hyperparam_interface.alpha_u);
            this.M = M;
            this.S = S;
            this.u_hyperparam_interface = u_hyperparam_interface;
            this.determine_u_hyperparams = MD_Determine_u_Hyperparameters(data_interface, u_hyperparam_interface);

            if this.u_hyperparam_interface.beta_u == 0.0
                this.determine_u_hyperparams.Determine_beta_u();
            end
            this.Set_beta_u(this.u_hyperparam_interface.beta_u);

            m = size(this.S, 1);
            u_vec = zeros(m, 1);
            if this.u_hyperparam_interface.gsvd_num_sing_vals == 0
                this.determine_u_hyperparams.Determine_GSVD_Hyperparameters();
            end
            this.Compute_E_u_Inverse_GSVD(this.u_hyperparam_interface.gsvd_num_sing_vals, this.u_hyperparam_interface.gsvd_oversampling, this.u_hyperparam_interface.gsvd_num_subspace_iter, u_vec);

            if ~this.u_hyperparam_interface.is_transient
                if this.u_hyperparam_interface.alpha_u == 0.0
                    this.determine_u_hyperparams.Determine_alpha_u(this);
                end
                this.Set_alpha_u(this.u_hyperparam_interface.alpha_u);
            end
        end

    end

end
