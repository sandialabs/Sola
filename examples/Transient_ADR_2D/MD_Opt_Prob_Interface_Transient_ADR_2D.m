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

            for k = 1:size(u_in, 2)
                u_in(:, k) = this.Reorder_u(u_in(:, k));
            end

            num_vecs = size(u_in, 2);
            full_dim = length(this.basis1.ybar);
            r1 = this.basis1.r;
            r2 = this.basis2.r;
            n_t = this.opt.con.n_t;
            u_in_1_tmp = zeros(r1 * n_t, num_vecs);
            u_in_2_tmp = zeros(r2 * n_t, num_vecs);
            u_in_1 = u_in(1:full_dim * n_t, :);
            u_in_2 = u_in((full_dim * n_t + 1):end, :);
            for k = 1:this.opt.con.n_t
                tmpk = u_in_1((1 + (k - 1) * full_dim):(k * full_dim), :);
                u_in_1_tmp((1 + (k - 1) * r1):(k * r1), :) = this.basis1.Compress(tmpk);

                tmpk = u_in_2((1 + (k - 1) * full_dim):(k * full_dim), :);
                u_in_2_tmp((1 + (k - 1) * r2):(k * r2), :) = this.basis2.Compress(tmpk);
            end

            u_in_tmp = [u_in_1_tmp; u_in_2_tmp];
            tmp = this.opt.con.c_u_Transpose_Inverse_Apply(u_in_tmp, this.u_current, z);
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
            u = this.Reorder_u(u);
            [~, grad_u, ~] = this.obj_hifi.J(u, z);
            grad_u = this.Reorder_u_Inverse(grad_u);
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            for k = 1:size(u_in, 2)
                u_in(:, k) = this.Reorder_u(u_in(:, k));
            end
            u = this.Reorder_u(u);
            u_out = this.obj_hifi.J_uu_Apply(u_in, u, z);
            for k = 1:size(u_in, 2)
                u_out(:, k) = this.Reorder_u_Inverse(u_out(:, k));
            end
        end

        function [u_reorder] = Reorder_u(this, u)
            u1 = u(1:this.obj_hifi.n_x * this.obj_hifi.n_t);
            u2 = u((this.obj_hifi.n_x * this.obj_hifi.n_t + 1):end);
            u1 = reshape(u1, this.obj_hifi.n_x, this.obj_hifi.n_t);
            u2 = reshape(u2, this.obj_hifi.n_x, this.obj_hifi.n_t);
            u_reorder = [u1; u2];
            u_reorder = u_reorder(:);
        end

        function [u] = Reorder_u_Inverse(this, u_reorder)
            u_reorder = reshape(u_reorder, 2 * this.obj_hifi.n_x, this.obj_hifi.n_t);
            u1 = u_reorder(1:this.obj_hifi.n_x, :);
            u2 = u_reorder((this.obj_hifi.n_x + 1):end, :);
            u = [u1(:); u2(:)];
        end

    end

end
