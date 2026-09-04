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
    %   X_s = Y_s R'
    %
    % and
    %
    %   Y_s ~ MN_{m,k}(0, W_u^{-1}, I_k),
    %   R R' approx Sigma_beta.
    %
    % The lazy_Y object provides consistent repeated applications of
    %
    %   gamma -> Y_s gamma,
    %   u     -> Y_s' u.
    %
    % Therefore,
    %
    %   X_s beta = Y_s (R' beta),
    %   X_s' u   = R (Y_s' u).

    properties
        R
        lazy_Y
    end

    methods

        function this = MD_Breve_Beta_Sampler(R, u_prior_interface, output_dim, tol)
            arguments
                R {mustBeNumeric}
                u_prior_interface MD_u_Prior_Interface
                output_dim (1, 1) {mustBeNumeric}
                tol (1, 1) {mustBeNumeric} = 1e-10
            end

            this.R = R;

            right_dim = size(R, 2);

            this.lazy_Y = MD_Lazy_Matrix_Normal_Operator( ...
                u_prior_interface, ...
                right_dim, ...
                output_dim, ...
                tol);
        end

        function u_out = Eval(this, beta)
            beta = beta(:);

            if length(beta) ~= size(this.R, 1)
                error('MD_Breve_Beta_Sampler::Eval beta has wrong dimension.');
            end

            gamma = this.R' * beta;
            u_out = this.lazy_Y.Forward_Apply(gamma);
        end

        function u_out = Apply_Jacobian(this, beta_in)
            beta_in = beta_in(:);

            if length(beta_in) ~= size(this.R, 1)
                error('MD_Breve_Beta_Sampler::Apply_Jacobian input has wrong dimension.');
            end

            gamma = this.R' * beta_in;
            u_out = this.lazy_Y.Forward_Apply(gamma);
        end

        function beta_out = Apply_Jacobian_Transpose(this, u_in)
            u_in = u_in(:);

            gamma_out = this.lazy_Y.Adjoint_Apply(u_in);
            beta_out = this.R * gamma_out;
        end

    end

end