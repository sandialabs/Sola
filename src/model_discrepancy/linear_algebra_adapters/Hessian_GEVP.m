%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Hessian_GEVP < Randomized_GEVP

    properties
        z
        opt_prob_interface
        z_prior_interface
        normalization_coeff
    end

    %% Implementation of base class virtual functions
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

        function [samples] = Generate_Random_Samples(this, num_samples)
            samples = this.z_prior_interface.Sample_with_Covariance_W_z_Inverse(num_samples);
            if size(samples, 1) > 0
                samples = sqrt(this.normalization_coeff) * samples;
            else
                samples = randn(length(this.z), num_samples);
            end
        end

    end

    %% Constructor and helper functions
    methods

        function this = Hessian_GEVP(opt_prob_interface, z_prior_interface, z)
            arguments
                opt_prob_interface MD_Opt_Prob_Interface
                z_prior_interface MD_z_Prior_Interface
                z (:, 1) double
            end
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
