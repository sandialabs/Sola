classdef MD_Opt_Prob_Interface_Transient_ADR_2D < MD_Opt_Prob_Interface

    properties
        opt
        obj_hifi
        basis1
        basis2
        z_current
        u_current
        hessian_data
        m
    end

    methods

        function this = MD_Opt_Prob_Interface_Transient_ADR_2D(opt, obj_hifi, basis1, basis2, z_lofi)
            this@MD_Opt_Prob_Interface();
            this.opt = opt;
            this.obj_hifi = obj_hifi;
            this.basis1 = basis1;
            this.basis2 = basis2;
            this.z_current = z_lofi;
            [~, ~, this.hessian_data] = this.opt.Jhat(this.z_current);
            this.m = (length(this.hessian_data) - length(this.z_current)) / 2;
            this.u_current = this.hessian_data(1:this.m);
        end

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.opt.Jhat(this.z_current);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m);
            end

            num_vecs = size(u_in, 2);
            u_in_red = zeros(this.m, num_vecs);
            full_dim = length(this.basis1.ybar);
            r1 = this.basis1.r;
            r2 = this.basis2.r;
            n_t = this.opt.con.n_t;
            for j = 1:num_vecs
                u_tmp = reshape(u_in(:, j), [], n_t);
                u_in_1 = zeros(r1, n_t);
                u_in_2 = zeros(r2, n_t);
                for k = 1:n_t
                    u_in_1(:, k) = this.basis1.V' * u_tmp(1:full_dim, k);
                    u_in_2(:, k) = this.basis2.V' * u_tmp((full_dim + 1):end, k);
                end
                tmp = [u_in_1; u_in_2];
                u_in_red(:, j) = tmp(:);
            end

            tmp = this.opt.con.c_u_Transpose_Inverse_Apply(u_in_red, this.u_current, z);
            z_out = -this.opt.con.c_z_Transpose_Apply(tmp, this.u_current, z);
        end

        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.opt.Jhat(this.z_current);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m);
            end
            z_out = this.opt.Jhat_hessVec(this.hessian_data, z_in);
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            [~, grad_u, ~] = this.obj_hifi.J(u, z);
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = this.obj_hifi.J_uu_Apply(u_in, u, z);
        end

    end

end
