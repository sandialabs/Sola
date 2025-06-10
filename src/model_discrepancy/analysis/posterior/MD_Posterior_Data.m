classdef MD_Posterior_Data < handle

    properties
        alpha_d
        N
        Z
        D
        Mz_Z
        Mz_z_opt
        Wz_inv_Mz_Z
        Wz_inv_Mz_z_opt
        Mz_Wz_inv_Mz_Z
        Mz_Wz_inv_Mz_z_opt
        Mz_Zc
        Zc
        Wz_inv_Mz_Zc
        Mz_Wz_inv_Mz_Zc
        Zc_Mz_Wz_inv_Mz_Zc
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
        Mz_z_breve
    end

    methods

        function this = MD_Posterior_Data()

        end

        function [] = Compute_Posterior_Data(this, data_interface, u_prior_interface, z_prior_interface, alpha_d_in, num_samples)
            this.alpha_d = alpha_d_in;
            this.num_samples = num_samples;
            this.z_opt = data_interface.z_opt;
            this.Z = data_interface.Z;
            this.D = data_interface.D;
            this.N = size(this.D, 2);

            this.Mz_Z = z_prior_interface.Apply_M_z(this.Z);
            this.Mz_z_opt = z_prior_interface.Apply_M_z(data_interface.z_opt);
            this.Wz_inv_Mz_Z = z_prior_interface.Apply_W_z_Inverse(this.Mz_Z);
            this.Wz_inv_Mz_z_opt = z_prior_interface.Apply_W_z_Inverse(this.Mz_z_opt);

            this.Mz_Wz_inv_Mz_Z = z_prior_interface.Apply_M_z(this.Wz_inv_Mz_Z);
            this.Mz_Wz_inv_Mz_z_opt = z_prior_interface.Apply_M_z(this.Wz_inv_Mz_z_opt);

            this.Mz_Zc = this.Mz_Z(:, 2:end) - this.Mz_z_opt;
            this.Wz_inv_Mz_Zc = this.Wz_inv_Mz_Z(:, 2:end) - this.Wz_inv_Mz_z_opt;
            this.Mz_Wz_inv_Mz_Zc = this.Mz_Wz_inv_Mz_Z(:, 2:end) - this.Mz_Wz_inv_Mz_z_opt;
            this.Zc_Mz_Wz_inv_Mz_Zc = this.Mz_Zc' * this.Wz_inv_Mz_Zc;
            this.Zc = this.Z(:, 2:end) - data_interface.z_opt;

            this.G = (1 + this.Wz_inv_Mz_z_opt' * this.Mz_z_opt) - this.Mz_Z' * this.Wz_inv_Mz_z_opt - this.Wz_inv_Mz_z_opt' * this.Mz_Z + this.Mz_Z' * this.Wz_inv_Mz_Z;
            if size(this.G, 1) == 1
                this.g_vecs = 1;
                this.Mu = 1;
            else
                [this.g_vecs, this.Mu] = eigs(this.G, this.N);
            end
            this.g_vecs = this.g_vecs .* (ones(this.N, 1) * sign(this.g_vecs(1, :)));

            M_u_D = u_prior_interface.Apply_M_u(this.D);
            this.u_ell = u_prior_interface.Apply_W_u_Inverse(M_u_D);
            this.u_i_ell = cell(this.N, 1);
            M_u_u_ell = u_prior_interface.Apply_M_u(this.u_ell);
            for i = 1:this.N
                this.u_i_ell{i} = (1 / this.alpha_d) * u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(M_u_u_ell, this.Mu(i, i) / this.alpha_d);
            end

            this.a_ell = zeros(this.N, 1);
            this.b_i_ell = zeros(this.N, this.N);
            for ell = 1:this.N
                this.a_ell(ell) = 1 - this.Wz_inv_Mz_z_opt' * (this.Mz_Z(:, ell) - this.Mz_z_opt);
                for i = 1:this.N
                    this.b_i_ell(i, ell) = (this.Mz_Z * this.g_vecs(:, i))' * (this.Wz_inv_Mz_Z(:, ell) - this.Wz_inv_Mz_z_opt) + sum(this.g_vecs(:, i)) * this.a_ell(ell);
                end
            end

            if this.num_samples > 0

                this.ui_hat = cell(this.N, 1);
                for i = 1:this.N
                    this.ui_hat{i} = (1 / sqrt(this.alpha_d)) * u_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this.num_samples, this.Mu(i, i) / this.alpha_d);
                end

                this.u_breve = u_prior_interface.Sample_with_Covariance_W_u_Inverse(this.num_samples);

                this.z_breve = z_prior_interface.Sample_with_Covariance_W_z_Inverse(this.num_samples);
                this.Mz_z_breve = z_prior_interface.Apply_M_z(this.z_breve);

            end

        end

    end

end
