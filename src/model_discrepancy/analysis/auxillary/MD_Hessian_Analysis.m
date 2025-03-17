classdef MD_Hessian_Analysis < handle

    properties
        opt_prob_interface
        z_prior_interface
        z_current
        evecs
        evals
    end

    methods

        function this = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface)
            this.opt_prob_interface = opt_prob_interface;
            this.z_prior_interface = z_prior_interface;
        end

        function [] = Compute_Hessian_GEVP(this, z, num_evals, oversampling)
            gevp = Hessian_GEVP(this.opt_prob_interface, this.z_prior_interface, z);
            [this.evecs, this.evals] = gevp.Compute_Hessian_GEVP(num_evals, oversampling);
            this.z_current = z;
        end

        function [z_out] = Apply_RS_Hessian_Inverse(this, z_in, z)
            if ~isempty(this.evals)
                if norm(this.z_current - z) ~= 0
                    disp('The z input has changed. Need to recompute the GEVP.');
                end
                z_out = this.Apply_Projected_RS_Hessian_Inverse(z_in);
            else
                z_out = this.Apply_RS_Hessian_Inverse_CG(z_in, z);
            end
        end

        function [z_out] = Apply_Projected_RS_Hessian_Inverse(this, z_in)
            z_out = this.evecs * diag(1 ./ this.evals) * this.evecs' * z_in;
        end

        function [z_out] = Apply_RS_Hessian_Inverse_CG(this, z_in, z)
            z_out = 0 * z_in;
            for k = 1:size(z_in, 2)
                tol = 1.e-7;
                max_iter = length(z);
                [z_out(:, k), flag, relres, iter, resvec] = pcg(@(x)this.opt_prob_interface.Apply_RS_Hessian(x, z), z_in(:, k), tol, max_iter);
                if flag ~= 0
                    disp('CG did not converge');
                end
            end
        end

    end

end
