classdef MD_Analytic_Laplacian_u_Prior_Interface < MD_Scaled_u_Prior_Interface

    properties
        M
        hyperparams
        beta_u
        sing_vecs_output
        sing_vals
    end

    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [u_out] = Apply_W_u_Acute_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            K = (this.sing_vals.^2) ./ (1 + scalar * this.sing_vals.^2);
            u_out = this.sing_vecs_output * diag(K) * this.sing_vecs_output' * u_in;
        end

        function [u_out] = Apply_W_u_Acute_Inverse(this, u_in)
            u_out = this.sing_vecs_output * diag(this.sing_vals.^2) * this.sing_vecs_output' * u_in;
        end

        % Compute samples from a mean zero Gaussian with covariance W_u^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Inverse(this, num_samples)
            r = length(this.sing_vals);
            u_out = sqrt(this.alpha_u) * this.sing_vecs_output * diag(this.sing_vals) * randn(r, num_samples);
        end

        % Compute samples from a mean zero Gaussian with covariance (W_u+scalar*M_u)^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            K = (this.sing_vals.^2) ./ (1 + this.alpha_u * scalar * this.sing_vals.^2);
            r = length(this.sing_vals);
            u_out = sqrt(this.alpha_u) * this.sing_vecs_output * diag(sqrt(K)) * randn(r, num_samples);
        end

    end

    methods

        function [] = Generate_Spectral_Decomposition(this)
            
            nodes = this.hyperparams.Load_Node_Data();
            m = size(nodes,1);
            d = size(nodes,2);
            r = this.hyperparams.gsvd_num_sing_vals;

            if d == 1

                L = max(nodes(:,1)) - min(nodes(:,1));
                x = nodes(:,1) - min(nodes(:,1));
                e = 1 + this.beta_u * (pi/L)^2 * (0:(r-1)).^2;
                this.sing_vals = 1./e';
                this.sing_vecs_output = zeros(m,r);
                v = ones(m,1);
                this.sing_vecs_output(:,1) = v/sqrt(v'*this.Apply_M_u(v));
                for k = 2:r
                    v = cos((k-1)*(pi/L)*x);
                    this.sing_vecs_output(:,k) = v/sqrt(v'*this.Apply_M_u(v));
                end

            elseif d == 2

                rsq = floor(sqrt(r));
                r = zeros(2,1);
                L = zeros(2,1);
                for j = 1:2
                    L(j) = max(nodes(:,j)) - min(nodes(:,j));
                    tmp = (L(j)/pi)^2 * (1/this.beta_u) * (1/this.hyperparams.W_u_inv_spectral_gap - 1);
                    r(j) = round(sqrt(tmp));
                    r(j) = min(r(j),rsq);
                end

                I = (0:(r(1)-1))';
                evali = kron(I.^2,ones(r(2),1));
                I = (0:(r(2)-1))';
                evalj = kron(ones(r(1),1),I.^2);
                this.sing_vals = 1./(1 + this.beta_u * pi^2 * (evali/L(1)^2 + evalj/L(2)^2));

                x = nodes(:,1) - min(nodes(:,1));
                y = nodes(:,2) - min(nodes(:,2));
                v1 = cos(pi*x*(0:(r(1)-1))/L(1));
                v2 = cos(pi*y*(0:(r(2)-1))/L(2));
                this.sing_vecs_output = zeros(m,r(1)*r(2));
                for i = 1:r(1)
                    I = (1 + (i-1)*r(2)):(i*r(2));
                    this.sing_vecs_output(:,I) = v2 .* (v1(:,i)*ones(1,r(1)));
                end

                [~,I] = sort(this.sing_vals,'descend');
                this.sing_vals = this.sing_vals(I);
                this.sing_vecs_output = this.sing_vecs_output(:,I);
                normalize = sqrt(this.sing_vecs_output'*this.Apply_M_u(this.sing_vecs_output));
                for l = 1:r(1)*r(2)
                    this.sing_vecs_output(:,l) = this.sing_vecs_output(:,l)/normalize(l,l);
                end

            elseif d == 3
                
                disp('MD_Analytic_Laplacian_u_Prior_Interface is not supported in 3D')
            
            end

        end
      
        function [] = Set_beta_u(this,beta_u_new)
            this.beta_u = beta_u_new;
            this.Generate_Spectral_Decomposition();
        end

        function this = MD_Analytic_Laplacian_u_Prior_Interface(M,hyperparams)
            this@MD_Scaled_u_Prior_Interface(hyperparams.alpha_u)
            this.M = M;
            this.hyperparams = hyperparams;
            this.beta_u = 0.0;

            if this.hyperparams.beta_u == 0.0
                this.hyperparams.Determine_beta_u();
            end
            this.Set_beta_u(this.hyperparams.beta_u);

            if this.hyperparams.gsvd_num_sing_vals == 0
                this.hyperparams.Determine_GSVD_Hyperparameters();
            end

            this.Generate_Spectral_Decomposition();

            if this.hyperparams.alpha_u == 0.0
                this.hyperparams.Determine_alpha_u(this);
            end
            this.Set_alpha_u(this.hyperparams.alpha_u);
        end

    end

end
