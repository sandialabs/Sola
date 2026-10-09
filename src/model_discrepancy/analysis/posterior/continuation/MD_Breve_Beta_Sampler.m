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
        materialized_X
        input_dim
        output_dim
    end

    methods

        function this = MD_Breve_Beta_Sampler(hessian_analysis, z_prior_interface, post_data, u_prior_interface, output_dim, tol, materialize_sample)
            arguments
                hessian_analysis MD_Hessian_Analysis
                z_prior_interface MD_z_Prior_Interface
                post_data MD_Posterior_Data
                u_prior_interface MD_u_Prior_Interface
                output_dim (1, 1) {mustBeNumeric}
                tol (1, 1) {mustBeNumeric} = 1e-10
                materialize_sample (1, 1) logical = false
            end
            this.hessian_analysis = hessian_analysis;
            this.z_prior_interface = z_prior_interface;
            this.post_data = post_data;

            if isempty(hessian_analysis.evals)
                right_dim = length(post_data.Mz_z_opt);
            else
                right_dim = length(hessian_analysis.evals);
            end

            this.input_dim = right_dim;
            this.output_dim = output_dim;
            this.materialized_X = [];

            sigma_apply = @(beta_in) this.Apply_Sigma_Beta(beta_in);
            sigma_sample = @(num_samples) this.Sample_Sigma_Beta(num_samples);
            this.lazy_X = MD_Lazy_Matrix_Normal_Operator(u_prior_interface, sigma_apply, sigma_sample, right_dim, output_dim, tol);

            if materialize_sample
                this.Materialize();
            end
        end

        function u_out = Eval(this, beta)
            beta = beta(:);
            if length(beta) ~= this.input_dim
                error('MD_Breve_Beta_Sampler::Eval beta has wrong dimension.');
            end
            if this.Is_Materialized()
                u_out = this.materialized_X * beta;
            else
                u_out = this.lazy_X.Forward_Apply(beta);
            end
        end

        function u_out = Apply_Jacobian(this, beta_in)
            beta_in = beta_in(:);
            if length(beta_in) ~= this.input_dim
                error('MD_Breve_Beta_Sampler::Apply_Jacobian input has wrong dimension.');
            end
            if this.Is_Materialized()
                u_out = this.materialized_X * beta_in;
            else
                u_out = this.lazy_X.Forward_Apply(beta_in);
            end
        end

        function beta_out = Apply_Jacobian_Transpose(this, u_in)
            u_in = u_in(:);
            if this.Is_Materialized()
                beta_out = this.materialized_X' * u_in;
            else
                beta_out = this.lazy_X.Adjoint_Apply(u_in);
            end
        end

        function [] = Materialize(this)
            if this.Is_Materialized()
                return;
            end

            X = zeros(this.output_dim, this.input_dim);
            I = eye(this.input_dim);
            for j = 1:this.input_dim
                X(:, j) = this.lazy_X.Forward_Apply(I(:, j));
            end
            this.materialized_X = X;
        end

        function [tf] = Is_Materialized(this)
            tf = ~isempty(this.materialized_X);
        end

        function [X] = Get_Materialized_Matrix(this)
            if ~this.Is_Materialized()
                this.Materialize();
            end
            X = this.materialized_X;
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
