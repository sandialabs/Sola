classdef MD_Laplacian_u_Prior_Interface < MD_Scaled_u_Prior_Interface

    properties
        M
        S

        M_lumped
        M_u_sqrt
        R_M
        P_M

        beta_u
        E_u
        is_sparse
        R_E
        P_E

        u_hyperparam_interface
        determine_u_hyperparams
    end

    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
             u_out = this.M * u_in;
        end

        function [u_out] = Apply_W_u_Acute_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            u_out = 0 * u_in;
            for k = 1:size(u_in, 2)
                tol = 1.e-7;
                max_iter = size(u_in,1);
                [u_out(:, k), flag, relres, iter, resvec] = pcg(@(x)this.Apply_W_u_Acute_Plus_scalar_M_u(x, scalar), u_in(:, k), tol, max_iter, @(x)this.Apply_W_u_Acute_Plus_scalar_M_u_Inverse_Approx(x, scalar));
                if flag ~= 0
                    disp('CG did not converge');
                end
            end
        end

        function [u_out] = Apply_W_u_Acute_Inverse(this, u_in)
            tmp1 = this.Apply_E_u_Inverse_Transpose(u_in);
            tmp2 = this.Apply_M_u(tmp1);
            u_out = this.Apply_E_u_Inverse(tmp2);
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Inverse(this, num_samples)
            omega = randn(size(this.M, 1), num_samples);
            vec = this.M_u_sqrt.Matrix_Sqrt_Apply(omega);
            u_out = this.Apply_E_u_Inverse(vec);
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(this, num_samples, scalar)

            n_u = size(this.M,1);
            u_out = zeros(n_u,num_samples);

            A = this.E_u' * diag(1./this.M_lumped) * this.E_u + scalar * this.M;
            R = chol(A);
            u_trial = linsolve(R,randn(n_u,num_samples));

            tmp1 = this.Apply_W_u_Acute_Plus_scalar_M_u(u_trial,scalar) - A*u_trial;
            tmp2 = sum(u_trial' .* permute(tmp1, [2, 1, 3]), 2); % Forming u_trial' * A * u_trial
            r = exp(-0.5 * tmp2);
            u = rand(num_samples,1);
            I = find(u < r);

            if max(r) > 1
                disp('Error: max(r) = should not exceed 1')
            end

            samples_so_far = length(I);
            u_out(:,1:samples_so_far) = u_trial(:,I);

            accept_rate = 1-(num_samples-samples_so_far)/num_samples;
            if (accept_rate < 0.25) || (accept_rate > 0.9)
                disp(['Warning: Acceptance rate in Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse is ',num2str(accept_rate)])
            end

            samples_to_generate = round((num_samples-samples_so_far)/accept_rate);
            while samples_to_generate > 0
                u_trial = linsolve(R,randn(n_u,samples_to_generate));
                tmp1 = this.Apply_W_u_Acute_Plus_scalar_M_u(u_trial,scalar) - A*u_trial;
                tmp2 = sum(u_trial' .* permute(tmp1, [2, 1, 3]), 2); % Forming u_trial' * A * u_trial
                r = exp(-0.5 * tmp2);
                u = rand(samples_to_generate,1);

                I = find(u < r);
                samples_here = length(I);
                if samples_here > num_samples-samples_so_far
                    I = I(1:(num_samples-samples_so_far));
                    samples_here = num_samples-samples_so_far;
                end

                u_out(:,(samples_so_far+1):(samples_so_far+samples_here)) = u_trial(:,I);
                samples_so_far = samples_so_far+samples_here;
                samples_to_generate = num_samples - samples_so_far;
            end
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface)
            arguments
                S (:, :) {mustBeNumeric}
                M (:, :) {mustBeNumeric}
                data_interface MD_Data_Interface
                u_hyperparam_interface MD_u_Hyperparameter_Interface
            end
            this@MD_Scaled_u_Prior_Interface(u_hyperparam_interface.alpha_u);

            this.M = M;
            this.S = S;
            this.is_sparse = issparse(this.M);

            this.M_lumped = sum(this.M,1);
            this.M_u_sqrt = M_u_Sqrt(this);
            if this.is_sparse
                [this.R_M, flag, this.P_M] = chol(this.M);
                if flag ~= 0
                    disp('Error in Cholesky factorization of M');
                end
            else
                this.R_M = chol(this.M);
            end

            this.u_hyperparam_interface = u_hyperparam_interface;
            this.determine_u_hyperparams = MD_Determine_u_Hyperparameters(data_interface, u_hyperparam_interface);

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
            this.E_u = this.beta_u * this.S + this.M;

            if this.is_sparse
                [this.R_E, flag, this.P_E] = chol(this.E_u);
                if flag ~= 0
                    disp('Error in Cholesky factorization of E_u');
                end
            else
                this.R_E = chol(this.E_u);
            end
        end

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            if this.is_sparse
                tmp1 = this.R_E' \ (this.P_E' * u_in);
                tmp2 = this.R_E \ tmp1;
                u_out = this.P_E * tmp2;
            else
                tmp = this.R_E' \ u_in;
                u_out = this.R_E \ tmp;
            end
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            if this.is_sparse
                tmp1 = this.R_E' \ (this.P_E' * u_in);
                tmp2 = this.R_E \ tmp1;
                u_out = this.P_E * tmp2;
            else
                tmp = this.R_E' \ u_in;
                u_out = this.R_E \ tmp;
            end
        end

        function [u_out] = Apply_M_u_Inverse(this, u_in)
            if this.is_sparse
                tmp1 = this.R_M' \ (this.P_M' * u_in);
                tmp2 = this.R_M \ tmp1;
                u_out = this.P_M * tmp2;
            else
                tmp = this.R_M' \ u_in;
                u_out = this.R_M \ tmp;
            end
        end

        function [u_out] = Apply_W_u_Acute(this, u_in)
            u_out = this.E_u' * this.Apply_M_u_Inverse(this.E_u * u_in);
        end

        function [u_out] = Apply_W_u_Acute_Plus_scalar_M_u(this, u_in, scalar)
            u_out = this.Apply_W_u_Acute(u_in) + scalar * this.Apply_M_u(u_in);
        end

        function [u_out] = Apply_W_u_Acute_Plus_scalar_M_u_Inverse_Approx(this, u_in, scalar)
            A = this.E_u' * diag(1./this.M_lumped) * this.E_u + scalar * this.M;
            u_out = A \ u_in;
        end
    end

end
