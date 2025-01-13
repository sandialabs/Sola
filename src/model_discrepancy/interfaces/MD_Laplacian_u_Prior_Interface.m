classdef MD_Laplacian_u_Prior_Interface < MD_Elliptic_u_Prior_Interface

    properties
        M
        S
        hyperparams
        E_u
        beta_u
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = this.E_u \ u_in;
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = this.E_u' \ u_in;
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [] = Set_beta_u(this,beta_u_new)
            this.beta_u = beta_u_new;
            this.E_u = this.beta_u * this.S + this.M;
        end

        function this = MD_Laplacian_u_Prior_Interface(S,M,hyperparams)
            this@MD_Elliptic_u_Prior_Interface(hyperparams.alpha_u)
            this.M = M;
            this.S = S;
            this.hyperparams = hyperparams;

            if this.hyperparams.beta_u == 0.0
                this.hyperparams.Determine_beta_u();
            end
            this.Set_beta_u(this.hyperparams.beta_u);

            m = size(this.S,1);
            u_vec = zeros(m,1);
            if this.hyperparams.gsvd_num_sing_vals == 0
                this.hyperparams.Determine_GSVD_Hyperparameters(m);
            end
            this.Compute_E_u_Inverse_GSVD(this.hyperparams.gsvd_num_sing_vals, this.hyperparams.gsvd_oversampling, this.hyperparams.gsvd_num_subspace_iter, u_vec);

            if this.hyperparams.alpha_u == 0.0
                this.hyperparams.Determine_alpha_u(this);
            end
            this.Set_alpha_u(this.hyperparams.alpha_u);
        end

    end

end
