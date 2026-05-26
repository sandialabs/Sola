%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Linear_OED_D_Opt < handle

    % We assume a linear Bayesian inverse problem with a mean zero Gaussian
    % noise, a Gaussian prior, and a linear observation operator
    %
    % Requires cons.c_u_Transpose_Inverse_Apply is the adjoint with respect
    % to the appropriate weighted inner product.
    %
    % Assumes a cardinality constraint on the design and uses a greedy algorithm to
    % determine an design.

    properties
        likelihood
        prior
        con
        num_sensors

        u_dim
        z_dim
        d_dim
        sigma

        Fm_cov
        G_noise
    end

    methods (Access = public)

        function this = Linear_OED_D_Opt(likelihood, prior, con, num_sensors)
            this.likelihood = likelihood;
            this.prior = prior;
            this.con = con;
            this.num_sensors = num_sensors;
            z = prior.Get_Prior_Mean();
            this.z_dim = length(z);
            u = con.c_z_Apply(z);
            this.u_dim = length(u);
            d = likelihood.Observation_Operator_Apply(u);
            this.d_dim = length(d);

            % Construct the operator F. This assumes that the data dimension
            % is very small in comparison to the parameter dimension
            Fm_cov = zeros(this.d_dim, this.d_dim);
            for i = 1:this.d_dim
                temp = zeros(this.d_dim, 1);
                temp(i) = 1.0;
                temp = this.likelihood.Observation_Operator_Transpose_Apply(temp);
                % Compute adjoint of solution map
                temp = this.con.c_u_Transpose_Inverse_Apply(temp);
                temp = -this.con.c_z_Transpose_Apply(temp);
                temp = this.prior.Mass_Matrix_Inverse_Apply(temp);
                temp = this.prior.Prior_Covariance_Apply(temp);
                temp = this.con.State_Solve(temp);
                Fm_cov(i, :) = this.likelihood.Observation_Operator_Apply(temp);
            end
            this.Fm_cov = Fm_cov;

            % Construct the noise covariance in a brute force fashion
            noise_prec = zeros(this.d_dim, this.d_dim);
            for i = 1:this.d_dim
                x = zeros(this.d_dim, 1);
                x(i) = 1.0;
                noise_prec(:, i) = this.likelihood.Noise_Precision_Apply(x);
            end
            this.G_noise = inv(noise_prec);
        end

        function [val] = OED_Objective(this, w)
            temp = det((this.Fm_cov(w, w)) + this.G_noise(w, w)) / det(this.G_noise(w, w));
            val = 0.5 * log(temp);
        end

        function [w] = Optimize_Design(this)
            w = Lazy_Greedy_Solve_Cardinality_Cons(@(S) this.OED_Objective(S), ...
                                                   this.d_dim, ...
                                                   this.num_sensors);
        end

    end

end
