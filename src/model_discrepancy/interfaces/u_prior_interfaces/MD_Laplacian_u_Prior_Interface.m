classdef MD_Laplacian_u_Prior_Interface < MD_Scaled_u_Prior_Interface

    properties
        M
        S

        use_lumped_mass
        M_lumped
        M_sqrt

        beta_u
        E_u
        is_sparse
        R_E
        P_E
        R_M
        P_M

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
                [u_out(:, k), flag, relres, iter, resvec] = pcg(@(x)this.Apply_W_u_Acute_Plus_scalar_M_u(x, scalar), u_in(:, k), tol, max_iter, @(x)this.Apply_W_u_Acute_Inverse(x));
                if flag ~= 0
                    disp('CG did not converge');
                end
            end
        end

        function [u_out] = Apply_W_u_Acute_Inverse(this, u_in)
            tmp1 = this.Apply_E_u_Inverse_Transpose(u_in);
            if this.use_lumped_mass
                tmp2 = diag(this.M_lumped) * tmp1;
            else
                tmp2 = this.Apply_M_u(tmp1);
            end
            u_out = this.Apply_E_u_Inverse(tmp2);
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Inverse(this, num_samples)
            omega = randn(size(this.M, 1), num_samples);
            if this.use_lumped_mass
                vec = diag(sqrt(this.M_lumped)) * omega;
            else
                vec = this.M_sqrt.Matrix_Sqrt_Apply(omega);
            end
            u_out = this.Apply_E_u_Inverse(vec);
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(this, num_samples, scalar)

            u_out = zeros(size(this.M,1),num_samples);

            u_trial = this.Sample_with_Covariance_W_u_Acute_Inverse(num_samples);
            tmp1 = this.Apply_M_u(u_trial);
            tmp2 = sum(u_trial' .* permute(tmp1, [2, 1, 3]), 2); % Forming u_trial' * M * u_trial
            r = exp(-0.5 * scalar * tmp2);
            u = rand(num_samples,1);
            I = find(u < r);

            samples_so_far = length(I);
            u_out(:,1:samples_so_far) = u_trial(:,I);

            accept_rate = 1-(num_samples-samples_so_far)/num_samples;
            if accept_rate < 0.8
                disp(['Warning: Acceptance rate in Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse is ',num2str(accept_rate)])
            end

            samples_to_generate = round((num_samples-samples_so_far)/accept_rate);
            while samples_to_generate > 0
                u_trial = this.Sample_with_Covariance_W_u_Acute_Inverse(samples_to_generate);
                tmp1 = this.Apply_M_u(u_trial);
                tmp2 = sum(u_trial' .* permute(tmp1, [2, 1, 3]), 2); % Forming u_trial' * M * u_trial
                r = exp(-0.5 * scalar * tmp2);
                u = rand(samples_to_generate,1);

                I = find(u < r);
                samples_here = length(I);
                if samples_here > samples_to_generate
                    I = I(1:samples_to_generate);
                    samples_here = samples_to_generate;
                end

                u_out(:,(samples_so_far+1):(samples_so_far+samples_here)) = u_trial(:,I);
                samples_so_far = samples_so_far+samples_here;
                samples_to_generate = num_samples - samples_so_far;
            end
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface, use_lumped_mass)
            arguments
                S (:, :) {mustBeNumeric}
                M (:, :) {mustBeNumeric}
                data_interface MD_Data_Interface
                u_hyperparam_interface MD_u_Hyperparameter_Interface
                use_lumped_mass {boolean} = false
            end
            this@MD_Scaled_u_Prior_Interface(u_hyperparam_interface.alpha_u);

            this.M = M;
            this.S = S;
            this.use_lumped_mass = use_lumped_mass;

            if use_lumped_mass
                this.M_lumped = sum(this.M,1);
            else
                this.M_sqrt = M_u_Sqrt(this);
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
            this.is_sparse = issparse(this.E_u);

            if this.is_sparse
            
                [this.R_E, flag, this.P_E] = chol(this.E_u);
                if flag ~= 0
                    disp('Error in Cholesky factorization of E_u');
                end

                if ~this.use_lumped_mass
                    [this.R_M, flag, this.P_M] = chol(this.M);
                    if flag ~= 0
                        disp('Error in Cholesky factorization of M');
                    end
                end
            
            else
            
                this.R_E = chol(this.E_u);
                if ~this.use_lumped_mass
                    this.R_M = chol(this.M);
                end
            
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
            if this.use_lumped_mass
                u_out = diag(1./this.M_lumped) * u_in;
            else

                if this.is_sparse
                    tmp1 = this.R_M' \ (this.P_M' * u_in);
                    tmp2 = this.R_M \ tmp1;
                    u_out = this.P_M * tmp2;
                else
                    tmp = this.R_M' \ u_in;
                    u_out = this.R_M \ tmp;
                end

            end
        end

        function [u_out] = Apply_W_u_Acute(this, u_in)
            u_out = this.E_u' * this.Apply_M_u_Inverse(this.E_u * u_in);
        end

        function [u_out] = Apply_W_u_Acute_Plus_scalar_M_u(this, u_in, scalar)
            u_out = this.Apply_W_u_Acute(u_in) + scalar * this.Apply_M_u(u_in);
        end

    end

end
