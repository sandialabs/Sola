classdef Misfit_Hessian_GEVP < Randomized_GEVP

    properties
        forward_operator_sing_vecs_input
        forward_operator_sing_vecs_output
        forward_operator_sing_vals
        w
        inf_dim_prior
        sigma_sq
    end

    methods (Access = public)

        function [vec_out] = Apply_Operator(this, vec_in)
            tmp1 = this.inf_dim_prior.Mass_Matrix_Inverse_Apply(vec_in);
            tmp2 = this.inf_dim_prior.Laplacian_Like_Transpose_Apply(tmp1);
            tmp3 = this.forward_operator_sing_vecs_input' * tmp2;
            tmp4 = diag(this.forward_operator_sing_vals) * tmp3;
            tmp5 = this.forward_operator_sing_vecs_output * tmp4;
            tmp6 = diag(this.w) * tmp5;
            tmp7 = this.forward_operator_sing_vecs_output' * tmp6;
            tmp8 = diag(this.forward_operator_sing_vals) * tmp7;
            tmp9 = this.forward_operator_sing_vecs_input * tmp8;
            tmp10 = this.inf_dim_prior.Laplacian_Like_Apply(tmp9);
            vec_out = (1 / this.sigma_sq) * this.inf_dim_prior.Mass_Matrix_Inverse_Apply(tmp10);
        end

        function [vec_out] = Apply_Weighting_Operator(this, vec_in)
            vec_out = this.inf_dim_prior.Mass_Matrix_Inverse_Apply(vec_in);
        end

        function [vec_out] = Apply_Weighting_Operator_Inverse(this, vec_in)
            vec_out = this.inf_dim_prior.Mass_Matrix_Apply(vec_in);
        end

        function [samples] = Generate_Random_Samples(this, num_samples)
            samples = randn(size(this.forward_operator_sing_vecs_input, 1), num_samples);
        end

    end

    methods

        function this = Misfit_Hessian_GEVP(forward_operator_sing_vecs_input, forward_operator_sing_vecs_output, forward_operator_sing_vals, w, inf_dim_prior, sigma_sq)
            this@Randomized_GEVP(forward_operator_sing_vecs_input(:, 1));
            this.forward_operator_sing_vecs_input = forward_operator_sing_vecs_input;
            this.forward_operator_sing_vecs_output = forward_operator_sing_vecs_output;
            this.forward_operator_sing_vals = forward_operator_sing_vals;
            this.w = w;
            this.inf_dim_prior = inf_dim_prior;
            this.sigma_sq = sigma_sq;
        end

    end

end
