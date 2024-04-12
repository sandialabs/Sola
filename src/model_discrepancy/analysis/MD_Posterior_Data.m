classdef MD_Posterior_Data < handle

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
        u_ell
        u_i_ell
        a_ell
        b_i_ell

        num_samples
        ui_hat
        u_breve
        z_breve
    end

    methods

        function this = MD_Posterior_Data()

        end

        function [] = Compute_Posterior_Data(this, data_interface, u_prior_interface, z_prior_interface, alpha_d_in, z_opt, num_samples)
            this.alpha_d = alpha_d_in;
            this.num_samples = num_samples;
            this.Z = data_interface.Z;
            this.D = data_interface.D;
            this.N = size(this.D, 2);

            this.W_z_inv_Z = z_prior_interface.Apply_W_z_Inverse(this.Z);
            this.W_z_inv_z_opt = z_prior_interface.Apply_W_z_Inverse(z_opt);
            this.G = (1 + this.W_z_inv_z_opt' * z_opt) - this.Z' * this.W_z_inv_z_opt - this.W_z_inv_z_opt' * this.Z + this.Z' * this.W_z_inv_Z;
            [this.g_vecs, this.Mu] = eig(this.G);

            M_u_Y = u_prior_interface.Apply_M_u(this.D);
            this.u_ell = u_prior_interface.Apply_W_u_Inverse(M_u_Y);
            this.u_i_ell = cell(this.N, 1);
            M_u_u_ell = u_prior_interface.Apply_M_u(this.u_ell);
            for i = 1:this.N
                this.u_i_ell{i} = (1 / this.alpha_d) * u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(M_u_u_ell, this.Mu(i, i) / this.alpha_d);
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
                for i = 1:this.N
                    this.ui_hat{i} = (1 / sqrt(this.alpha_d)) * u_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this.num_samples, this.Mu(i, i) / this.alpha_d);
                end

                this.u_breve = u_prior_interface.Sample_with_Covariance_W_u_Inverse(this.num_samples);

                this.z_breve = z_prior_interface.Sample_with_Covariance_W_z_Inverse(this.num_samples);

            end

        end

    end

end
