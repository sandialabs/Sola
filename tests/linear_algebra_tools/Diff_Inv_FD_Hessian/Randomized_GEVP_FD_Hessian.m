classdef Randomized_GEVP_FD_Hessian < Randomized_GEVP

    properties
        M
        opt
        grad_nom
        z_star
    end

    methods (Access = public)

        function [vec_out] = Apply_Operator(this, vec_in)

            h = 1.e-4;
            vec_out = 0.0 * vec_in;
            for k = 1:size(vec_in,2)
                zk = this.z_star + h * vec_in(:,k);
                [~,gradk] = this.opt.Jhat(zk);
                vec_out(:,k) = (gradk - this.grad_nom)/h;
            end

        end

        function [vec_out] = Apply_Weighting_Operator(this, vec_in)
            vec_out = this.M * vec_in;
        end

        function [vec_out] = Apply_Weighting_Operator_Inverse(this, vec_in)
            vec_out = linsolve(this.M, vec_in);
        end

        function [samples] = Generate_Random_Samples(this, num_samples)
            R = chol(this.M);
            tmp = randn(size(R, 1), num_samples);
            samples = linsolve(R, tmp);
        end

    end

    methods

        function this = Randomized_GEVP_FD_Hessian(opt,z_star)
            m = size(opt.obj.prior.L,1);
            vec = zeros(m, 1);
            this@Randomized_GEVP(vec);
            this.M = opt.obj.prior.con.M;
            this.opt = opt;
            [~,this.grad_nom] = opt.Jhat(z_star);
            this.z_star = z_star;
        end

    end

end
