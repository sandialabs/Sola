classdef MD_Bayes_Posterior_Data < handle

    properties
        alpha_d
        N
        Z
        D
        W_z_inv_Z
        W_z_inv_z_opt
        G
        g_vecs
        Mu
        state_grad
        u_ell
        u_i_ell
        a_ell
        b_i_ell

        num_samples
        ui_hat
        u_breve
        state_grad_W_u_inv_state_grad
        zbreve
    end

    methods

        function this = MD_Bayes_Posterior_Data()

        end

        function [] = Compute_Posterior_Data(this, opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, alpha_d_in, u_opt, z_opt, num_samples)
            this.alpha_d = alpha_d_in;
            this.num_samples = num_samples;
            this.Z = data_interface.Z;
            this.D = data_interface.D;
            this.N = size(this.D, 2);
            this.state_grad = opt_prob_interface.Misfit_Gradient(u_opt, z_opt);

            this.W_z_inv_Z = z_prior_interface.Apply_W_z_Inverse(this.Z);
            this.W_z_inv_z_opt = z_prior_interface.Apply_W_z_Inverse(z_opt);
            this.G = (1 + this.W_z_inv_z_opt' * z_opt) - this.Z' * this.W_z_inv_z_opt - this.W_z_inv_z_opt' * this.Z + this.Z' * this.W_z_inv_Z;
            [this.g_vecs, this.Mu] = eig(this.G);

            W_d_Y = u_prior_interface.Apply_W_d(this.D);
            this.u_ell = u_prior_interface.Apply_W_u_Inverse(W_d_Y);
            this.u_i_ell = cell(this.N, 1);
            W_d_u_ell = u_prior_interface.Apply_W_d(this.u_ell);
            for i = 1:this.N
                this.u_i_ell{i} = (1 / this.alpha_d) * u_prior_interface.Apply_W_u_Plus_scalar_W_d_Inverse(W_d_u_ell, this.Mu(i, i) / this.alpha_d);
            end

            this.a_ell = zeros(this.N, 1);
            this.b_i_ell = zeros(this.N, this.N);
            for ell = 1:this.N
                this.a_ell(ell) = 1 - this.W_z_inv_z_opt' * (this.Z(:, ell) - z_opt);
                for i = 1:this.N
                    this.b_i_ell(i, ell) = (this.Z * this.g_vecs(:, i))' * (this.W_z_inv_Z(:, ell) - this.W_z_inv_z_opt) + sum(this.g_vecs(:, i)) * this.a_ell(ell);
                end
            end

            if this.num_samples > 0

                this.ui_hat = cell(this.N, 1);
                m = size(this.u_ell, 1);
                for i = 1:this.N
                    Omega = randn(m, this.num_samples);
                    this.ui_hat{i} = (1 / sqrt(this.alpha_d)) * u_prior_interface.Apply_W_u_Plus_scalar_W_d_Inverse_Factor(Omega, this.Mu(i, i) / this.alpha_d);
                end

                Omega = randn(m, this.num_samples);
                this.u_breve = u_prior_interface.Apply_W_u_Inverse_Factor(Omega);

                this.state_grad_W_u_inv_state_grad = u_prior_interface.Apply_W_u_Inverse(this.state_grad)' * this.state_grad;
                n = length(z_opt);
                Omega = randn(n, this.num_samples);
                this.zbreve = z_prior_interface.Apply_W_z_Inverse_Factor(Omega);

            end

        end

    end

end
