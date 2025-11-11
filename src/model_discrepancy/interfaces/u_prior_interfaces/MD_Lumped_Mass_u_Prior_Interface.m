classdef MD_Lumped_Mass_u_Prior_Interface < MD_Scaled_u_Prior_Interface

    properties
        M_lumped_diag
        M
        S
        E_u
        R
        P
        W_u_acute

        u_hyperparam_interface
        determine_u_hyperparams
        is_sparse
        use_sampling_prec
        beta_u
        n_u
    end

    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [u_out] = Apply_W_u_Acute_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            A = this.W_u_acute + scalar * this.M;
            u_out = A \ u_in;
        end

        function [u_out] = Apply_W_u_Acute_Inverse(this, u_in)
            tmp1 = this.Apply_E_u_Inverse_Transpose(u_in);
            tmp2 = diag(this.M_lumped_diag)*tmp1;
            u_out = this.Apply_E_u_Inverse(tmp2);
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Inverse(this, num_samples)
            omega = randn(this.n_u,num_samples);
            vec = diag(sqrt(this.M_lumped_diag)) * omega;
            u_out = this.Apply_E_u_Inverse(vec);
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            A = this.W_u_acute + scalar * this.M;
            if this.use_sampling_prec
                if this.is_sparse
                    L = ichol(A);
                else
                    L = chol(A)';
                end
                W_u_Acute_Plus_scalar_M_u_sqrt = Sparse_Matrix_Sqrt(A, L);
            else
                W_u_Acute_Plus_scalar_M_u_sqrt = Sparse_Matrix_Sqrt(A);
            end
            omega = randn(this.n_u,num_samples);
            tmp = W_u_Acute_Plus_scalar_M_u_sqrt.Matrix_Sqrt_Apply(omega);
            u_out = this.Apply_W_u_Acute_Plus_scalar_M_u_Inverse(tmp, scalar);
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Lumped_Mass_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface)
            arguments
                S (:, :) {mustBeNumeric}
                M (:, :) {mustBeNumeric}
                data_interface MD_Data_Interface
                u_hyperparam_interface MD_u_Hyperparameter_Interface
            end
            this@MD_Scaled_u_Prior_Interface(u_hyperparam_interface.alpha_u);

            this.M = M;
            this.S = S;
            this.u_hyperparam_interface = u_hyperparam_interface;
            this.determine_u_hyperparams = MD_Determine_u_Hyperparameters(data_interface, u_hyperparam_interface);
            this.is_sparse = issparse(this.M);
            this.n_u = size(this.M,1);
            this.use_sampling_prec = true;

            this.M_lumped_diag = this.M * ones(size(this.M,1),1);

            if this.u_hyperparam_interface.beta_u == 0.0
                this.determine_u_hyperparams.Determine_beta_u();
            end
            this.Set_beta_u(this.u_hyperparam_interface.beta_u);

            if ~this.u_hyperparam_interface.is_transient
                if this.u_hyperparam_interface.alpha_u == 0.0
                    this.determine_u_hyperparams.Determine_alpha_u(this);
                end
                this.Set_alpha_u(this.u_hyperparam_interface.alpha_u);
            end
        end

        function [] = Set_beta_u(this, beta_u_new)
            this.beta_u = beta_u_new;
            this.Assemble_Operators();
        end

        function [] = Assemble_Operators(this)
            this.E_u = this.beta_u * this.S + this.M;

            if this.is_sparse
                [this.R, flag, this.P] = chol(this.E_u);
                if flag ~= 0
                    disp('Error in Cholesky factorization of E_u');
                end
            else
                this.R = chol(this.E_u);
            end

            this.W_u_acute = this.E_u' * sparse(diag(1./this.M_lumped_diag)) * this.E_u;
        end

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

    end

end
