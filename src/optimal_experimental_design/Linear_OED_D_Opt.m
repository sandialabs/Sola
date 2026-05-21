%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Linear_OED_D_Opt < handle

    % We assume a linear Bayesian inverse problem with a mean zero Gaussian noise
    % model with covariance sigma^2*I, a Gaussian prior, and a linear observation operator
    %
    % Requires prior.Prior_Covariance_Factor_Apply to be the square root of the prior
    % covariance operator and that cons.c_u_Transpose_Inverse_Apply is the adjoint with
    % respect to the appropriate weighted inner product.
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

        F
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
            tmp = this.likelihood.Noise_Precision_Apply(ones(this.d_dim, 1));
            if max(tmp) - min(tmp) > 0
                disp('Error: the noise covariance must be a scalar multiple of the identity');
            end
            this.sigma = sqrt(1 / mean(tmp));

            % Construct the operator F. This assumes that the data dimension
            % is very small in comparison to the parameter dimension
            F = zeros(this.d_dim, this.z_dim);
            for i = 1:this.d_dim
                x = zeros(this.d_dim, 1);
                x(i) = 1;
                temp = this.likelihood.Observation_Operator_Transpose_Apply(x);
                temp = this.con.c_u_Transpose_Inverse_Apply(temp) * this.sigma;
                f = this.prior.Prior_Covariance_Factor_Apply(temp);
                F(i, :) = f;
            end
            this.F = F;
        end

        function [w, eig] = Optimize_Design(this)
            w = [];
            M = this.con.M;
            C = RankOneUpdatesMatrix(M);
            marginal_gains = RedBlackTree();
            eig = 0;

            for v = 1:this.d_dim
                marginal_gains.Insert(inf, ...
                                      struct('sensor', v, 'time', 0));
            end

            for i = 1:this.num_sensors
                while true
                    [gain, data] = marginal_gains.PopMax();
                    v = data.sensor;
                    t = data.time;
                    f = this.F(v, :)';
                    if t == i
                        w(i) = v;
                        C.Add_Update(f);
                        eig = eig + gain;
                        break
                    else
                        g = C.Inverse_Apply(f);
                        gain = log(1 + g' * M * f);
                        marginal_gains.Insert(gain, ...
                                              struct('sensor', v, 'time', i));
                    end
                end
            end
            eig = eig / 2;
        end

    end

end
