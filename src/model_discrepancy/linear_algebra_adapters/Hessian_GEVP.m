classdef Hessian_GEVP < Randomized_GEVP

    properties
        z
        opt_prob_interface
        z_prior_interface
        normalization_coeff
    end

    methods (Access = public)

        function [vec_out] = Apply_Operator(this, vec_in)
            vec_out = this.opt_prob_interface.Apply_RS_Hessian(vec_in, this.z);
        end

        function [vec_out] = Apply_Weighting_Operator(this, vec_in)
            vec_out = (1 / this.normalization_coeff) * this.z_prior_interface.Apply_W_z(vec_in);
        end

        function [vec_out] = Apply_Weighting_Operator_Inverse(this, vec_in)
            vec_out = this.normalization_coeff * this.z_prior_interface.Apply_W_z_Inverse(vec_in);
        end

        function [vec_out] = Apply_Weighting_Operator_Preconditioner_Factor(this, vec_in)
            vec_out = sqrt(this.normalization_coeff) * this.z_prior_interface.Apply_W_z_Inverse_Factor(vec_in);
        end

    end

    methods

        function this = Hessian_GEVP(opt_prob_interface, z_prior_interface, z)
            this@Randomized_GEVP(z);
            this.z = z;
            this.opt_prob_interface = opt_prob_interface;
            this.z_prior_interface = z_prior_interface;
            this.normalization_coeff = z_prior_interface.Apply_W_z(z)' * z;
        end

        function [evecs, evals] = Compute_Hessian_GEVP(this, num_evals, oversampling)
            [evecs, evals] = this.Compute_GEVP(num_evals, oversampling);
        end

    end

end
