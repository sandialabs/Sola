classdef MD_Transient_Vector_z_Prior_Interface < MD_Scaled_z_Prior_Interface

    properties
        beta_t
        S
        M
        n_t
        num_controls
        hyperparams
        E_t
        V
        Lambda
    end

    methods (Access = public)

        function [z_out] = Apply_M_z(this, z_in)
            z_out = kron(this.M,eye(this.num_controls)) * z_in;
        end

        function [z_out] = Apply_W_z_Acute_Inverse(this, z_in)
            z_out = linsolve(kron(this.E_t,eye(this.num_controls)),z_in);
        end

        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = linsolve(kron(this.M,eye(this.num_controls)),z_in);
        end

        function [z_out] = Sample_with_Covariance_W_z_Acute_Inverse(this, num_samples)
            z_out = zeros(this.n_t*this.num_controls,num_samples);
            for k = 1:this.num_controls
                I = k:this.num_controls:(this.n_t*this.num_controls);
                z_out(I,:) = this.V * sqrt(this.Lambda) * randn(this.n_t,num_samples);
            end
        end

        function [z_out] = Apply_W_z_Acute(this, z_in)
            z_out = kron(this.E_t ,eye(this.num_controls)) * z_in;
        end

        function [] = Set_beta_t(this, beta_t_new)
            this.beta_t = beta_t_new;
            this.E_t = this.beta_t * this.S + this.M;
            [this.V,this.Lambda] = eig(this.E_t);
        end

        function this = MD_Transient_Vector_z_Prior_Interface(S, M, num_controls, hyperparams)
            this@MD_Scaled_z_Prior_Interface(hyperparams.alpha_z);
            this.S = S;
            this.M = M;
            this.n_t = size(M,1);
            this.num_controls = num_controls;
            this.hyperparams = hyperparams;

            if hyperparams.beta_t == 0.0
                this.hyperparams.Determine_beta_t();
            end
            this.Set_beta_t(this.hyperparams.beta_t);

            if this.hyperparams.alpha_z == 0.0
                this.hyperparams.Determine_alpha_z(this);
            end
            this.Set_alpha_z(this.hyperparams.alpha_z);
        end

    end

end
