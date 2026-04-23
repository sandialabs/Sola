%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Bilaplacian_u_Prior_Interface < MD_Scaled_u_Prior_Interface

    properties
        M
        S
        E_u
        R_u
        R_mass

        M_lumped_diag
        W_u_acute_approx

        u_hyperparam_interface
        determine_u_hyperparams
        use_prec
        beta_u
        n_u

        W_u_Acute_Plus_scalar_M_u_Inverse_iters
        W_u_Acute_Plus_scalar_M_u_Sqrt_iters
    end

    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [u_out] = Apply_W_u_Acute_Plus_scalar_M_u_Inverse(this, u_in, scalar)

            if this.use_prec
                A_approx = this.W_u_acute_approx + scalar * this.M;
                R = chol(A_approx);
                Mfun = @(r) R \ (R' \ r);
            else
                Mfun = @(r) r;
            end

            Afun = @(r) this.Apply_W_u_Acute(r) + scalar * this.M * r;
            tol   = 1e-7;
            maxit = this.n_u;
            x0    = zeros(this.n_u,1);

            [u_out,flag,~,iter] = pcg(Afun, u_in, tol, maxit, Mfun, [], x0);
            if flag ~=0
                disp('Error in PCG')
            end
            this.W_u_Acute_Plus_scalar_M_u_Inverse_iters = [this.W_u_Acute_Plus_scalar_M_u_Inverse_iters ; iter];
        end

        function [u_out] = Apply_W_u_Acute_Inverse(this, u_in)
            tmp1 = this.Apply_E_u_Inverse_Transpose(u_in);
            tmp2 = this.Apply_M_u(tmp1);
            u_out = this.Apply_E_u_Inverse(tmp2);
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Inverse(this, num_samples)
            M_u_sqrt = Sparse_Matrix_Sqrt(this.M, this.R_mass');
            omega = randn(this.n_u, num_samples);
            vec = M_u_sqrt.Matrix_Sqrt_Apply(omega);
            u_out = this.Apply_E_u_Inverse(vec);
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            Afun = @(r) this.Apply_W_u_Acute(r) + scalar * this.M * r;
            if this.use_prec
                A_approx = this.W_u_acute_approx + scalar * this.M;
                L = ichol(A_approx);
            else
                L = eye(this.n_u);
            end
            W_u_Acute_Plus_scalar_M_u_sqrt = Sparse_Matrix_Sqrt(Afun, L);
            omega = randn(this.n_u, num_samples);
            [tmp,relres] = W_u_Acute_Plus_scalar_M_u_sqrt.Matrix_Sqrt_Apply(omega);
            this.W_u_Acute_Plus_scalar_M_u_Sqrt_iters = [this.W_u_Acute_Plus_scalar_M_u_Sqrt_iters ; length(relres{1})];
            u_out = this.Apply_W_u_Acute_Plus_scalar_M_u_Inverse(tmp, scalar);
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Bilaplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface)
            arguments
                S (:, :) {mustBeNumeric}
                M (:, :) {mustBeNumeric}
                data_interface MD_Data_Interface
                u_hyperparam_interface MD_u_Hyperparameter_Interface
            end
            this@MD_Scaled_u_Prior_Interface(u_hyperparam_interface.alpha_u);

            this.M = sparse(M);
            this.S = sparse(S);
            this.u_hyperparam_interface = u_hyperparam_interface;
            this.determine_u_hyperparams = MD_Determine_u_Hyperparameters(data_interface, u_hyperparam_interface);
            this.n_u = size(this.M, 1);
            this.use_prec = true;

            this.W_u_Acute_Plus_scalar_M_u_Inverse_iters = [];
            this.W_u_Acute_Plus_scalar_M_u_Sqrt_iters = [];

            this.M_lumped_diag = this.M * ones(this.n_u, 1);

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

            [this.R_u, flag] = chol(this.E_u);
            if flag ~= 0
                disp('Error in Cholesky factorization of E_u');
            end

            [this.R_mass, flag] = chol(this.M);
            if flag ~= 0
                disp('Error in Cholesky factorization of M_u');
            end

            this.W_u_acute_approx = this.E_u' * sparse(diag(1 ./ this.M_lumped_diag)) * this.E_u;

        end

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            tmp1 = this.R_u' \ u_in;
            u_out = this.R_u \ tmp1;
        end

        function [u_out] = Apply_M_u_Inverse(this, u_in)
            tmp1 = this.R_mass' \ u_in;
            u_out = this.R_mass \ tmp1;
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            tmp1 = this.R_u' \ u_in;
            u_out = this.R_u \ tmp1;
        end

        function [u_out] = Apply_W_u_Acute(this, u_in)
            tmp1 = this.E_u * u_in;
            tmp2 = this.Apply_M_u_Inverse(tmp1);
            u_out = this.E_u' * tmp2;
        end

    end

end
