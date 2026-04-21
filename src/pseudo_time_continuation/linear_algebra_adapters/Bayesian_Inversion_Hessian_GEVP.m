%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Bayesian_Inversion_Hessian_GEVP < Randomized_GEVP

    properties
        bayes_inv
        Gauss_Newton_Hess
        z
        theta
        u
        lambda
    end

    methods (Access = public)

        function [vec_out] = Apply_Operator(this, vec_in)
            w = this.bayes_inv.con.c_z_Apply(vec_in, this.u, this.z);
            mu = this.bayes_inv.con.c_u_Inverse_Apply(-w, this.u, this.z);
            yJ = this.bayes_inv.obj.J_uu_Apply(mu, this.u, this.z);
            if this.Gauss_Newton_Hess
                gamma = this.bayes_inv.con.c_u_Transpose_Inverse_Apply(-yJ, this.u, this.z);
                vec_out = this.bayes_inv.con.c_z_Transpose_Apply(gamma, this.u, this.z);
            else
                yc = this.bayes_inv.con.c_uu_Apply(mu, this.u, this.z, this.lambda) + this.bayes_inv.con.c_uz_Apply(vec_in, this.u, this.z, this.lambda);
                gamma = this.bayes_inv.con.c_u_Transpose_Inverse_Apply(-(yJ + yc), this.u, this.z);
                vec_out = this.bayes_inv.con.c_z_Transpose_Apply(gamma, this.u, this.z) + this.bayes_inv.con.c_zu_Apply(mu, this.u, this.z, this.lambda) + this.bayes_inv.con.c_zz_Apply(vec_in, this.u, this.z, this.lambda);
            end
        end

        function [vec_out] = Apply_Weighting_Operator(this, vec_in)
            vec_out = this.bayes_inv.prior.Prior_Precision_Apply(vec_in);
        end

        function [vec_out] = Apply_Weighting_Operator_Inverse(this, vec_in)
            vec_out = this.bayes_inv.prior.Prior_Covariance_Apply(vec_in);
        end

        function [samples] = Generate_Random_Samples(this, num_samples)
            z_in = randn(length(this.z), num_samples);
            samples = this.bayes_inv.prior.Prior_Covariance_Factor_Apply(z_in);
        end

    end

    methods

        function this = Bayesian_Inversion_Hessian_GEVP(z, theta, bayes_inv)
            this@Randomized_GEVP(z);
            this.z = z;
            this.theta = theta;
            this.bayes_inv = bayes_inv;
            this.Gauss_Newton_Hess = false;
        end

        function [evecs, evals, u, z, lambda, theta] = Compute_Hessian_GEVP(this, num_evals, oversampling)
            this.u = this.bayes_inv.con.State_Solve(this.z);
            [~, grad_u] = this.bayes_inv.obj.J(this.u, this.z);
            this.lambda = this.bayes_inv.con.c_u_Transpose_Inverse_Apply(-grad_u, this.u, this.z);
            [evecs, evals] = this.Compute_GEVP(num_evals, oversampling);
            u = this.u;
            z = this.z;
            lambda = this.lambda;
            theta = this.theta;
        end

    end

end
