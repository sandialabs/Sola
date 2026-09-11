%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Breve_Beta_Sampler < handle

    % Represents one posterior sample's beta-space breve discrepancy operator
    %
    %   beta -> X_s beta,
    %
    % where
    %
    %   X_s ~ MN_{m,r}(0, W_u^{-1}, Sigma_beta).
    %

    properties
        hessian_analysis
        z_prior_interface
        post_data
        lazy_X
    end

    methods

        function this = MD_Breve_Beta_Sampler(hessian_analysis, z_prior_interface, post_data, u_prior_interface, output_dim, tol)
            arguments
                hessian_analysis MD_Hessian_Analysis
                z_prior_interface MD_z_Prior_Interface
                post_data MD_Posterior_Data
                u_prior_interface MD_u_Prior_Interface
                output_dim (1, 1) {mustBeNumeric}
                tol (1, 1) {mustBeNumeric} = 1e-10
            end
            this.hessian_analysis = hessian_analysis;
            this.z_prior_interface = z_prior_interface;
            this.post_data = post_data;

            if isempty(hessian_analysis.evals)
                right_dim = length(post_data.Mz_z_opt);
            else
                right_dim = length(hessian_analysis.evals);
            end

            sigma_apply = @(beta_in) this.Apply_Sigma_Beta(beta_in);
            sigma_sample = @(num_samples) this.Sample_Sigma_Beta(num_samples);
            this.lazy_X = MD_Lazy_Matrix_Normal_Operator(u_prior_interface, sigma_apply, sigma_sample, right_dim, output_dim, tol);
        end

        function u_out = Eval(this, beta)
            beta = beta(:);
            if length(beta) ~= this.lazy_X.input_dim
                error('MD_Breve_Beta_Sampler::Eval beta has wrong dimension.');
            end
            u_out = this.lazy_X.Forward_Apply(beta);
        end

        function u_out = Apply_Jacobian(this, beta_in)
            beta_in = beta_in(:);
            if length(beta_in) ~= this.lazy_X.input_dim
                error('MD_Breve_Beta_Sampler::Apply_Jacobian input has wrong dimension.');
            end
            u_out = this.lazy_X.Forward_Apply(beta_in);
        end

        function beta_out = Apply_Jacobian_Transpose(this, u_in)
            u_in = u_in(:);
            beta_out = this.lazy_X.Adjoint_Apply(u_in);
        end

        function beta_out = Apply_Sigma_Beta(this, beta_in)
            beta_in = beta_in(:);
            z_in = this.hessian_analysis.Apply_V(beta_in);
            Mz_z_in = this.z_prior_interface.Apply_M_z(z_in);
            tmp_rhs = this.z_prior_interface.Apply_W_z_Inverse(Mz_z_in);

            if ~isempty(this.post_data.Zc_Mz_Wz_inv_Mz_Zc)
                coeff = linsolve(this.post_data.Zc_Mz_Wz_inv_Mz_Zc, this.post_data.Mz_Zc' * tmp_rhs);
                tmp_rhs = tmp_rhs - this.post_data.Wz_inv_Mz_Zc * coeff;
            end

            beta_out = this.hessian_analysis.Apply_V_Transpose(this.z_prior_interface.Apply_M_z(tmp_rhs));
        end

        function beta_samps = Sample_Sigma_Beta(this, num_samples)
            z_samps = this.z_prior_interface.Sample_with_Covariance_W_z_Inverse(num_samples);

            if ~isempty(this.post_data.Zc_Mz_Wz_inv_Mz_Zc)
                coeff = linsolve(this.post_data.Zc_Mz_Wz_inv_Mz_Zc, this.post_data.Mz_Zc' * z_samps);
                z_samps = z_samps - this.post_data.Wz_inv_Mz_Zc * coeff;
            end

            beta_samps = this.hessian_analysis.Apply_V_Transpose(this.z_prior_interface.Apply_M_z(z_samps));
        end

    end

end
