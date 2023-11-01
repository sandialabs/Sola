classdef Hessian_GEVP < Randomized_GEVP

    properties
        z_opt
        md_interface
        evals
        evecs
        is_computed
        normalization_coeff
    end

    methods (Access = public)

        function [vec_out] = Apply_Operator(this, vec_in)
            vec_out = this.md_interface.Apply_RS_Hessian(vec_in, this.z_opt);
        end

        function [vec_out] = Apply_Weighting_Operator(this, vec_in)
            vec_out = (1 / this.normalization_coeff) * this.md_interface.Apply_W_z(vec_in);
        end

        function [vec_out] = Apply_Weighting_Operator_Inverse(this, vec_in)
            vec_out = this.normalization_coeff * this.md_interface.Apply_W_z_Inverse(vec_in);
        end

        function [vec_out] = Apply_Weighting_Operator_Inverse_Factor(this, vec_in)
            vec_out = sqrt(this.normalization_coeff) * this.md_interface.Apply_W_z_Inverse_Factor(vec_in);
        end

    end

    methods

        function this = Hessian_GEVP(md_interface, z_opt)
            this@Randomized_GEVP(z_opt);
            this.z_opt = z_opt;
            this.md_interface = md_interface;
            this.is_computed = false;
            this.normalization_coeff = md_interface.Apply_W_z(z_opt)' * z_opt;
        end

        function [] = Compute_Hessian_GEVP(this, num_evals, oversampling)
            [this.evecs, this.evals] = this.Compute_GEVP(num_evals, oversampling);
            this.is_computed = true;
        end

        function [z_out] = Apply_Projected_RS_Hessian_Inverse(this, z_in)
            z_out = this.evecs * diag(1 ./ this.evals) * this.evecs' * z_in;
        end

    end

end
