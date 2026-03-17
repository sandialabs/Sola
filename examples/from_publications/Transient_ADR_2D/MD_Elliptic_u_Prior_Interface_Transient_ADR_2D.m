classdef MD_Elliptic_u_Prior_Interface_Transient_ADR_2D < MD_Elliptic_u_Prior_Interface

    properties
        M_u
        E_u
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            if 2 * size(u_in, 1) == size(this.E_u, 1)
                m = size(u_in, 1);
                u_out = this.E_u(1:m, 1:m) \ u_in;
            else
                u_out = this.E_u \ u_in;
            end
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            if 2 * size(u_in, 1) == size(this.E_u, 1)
                m = size(u_in, 1);
                u_out = this.E_u(1:m, 1:m)' \ u_in;
            else
                u_out = this.E_u' \ u_in;
            end
        end

        function [u_out] = Apply_M_u(this, u_in)
            if 2 * size(u_in, 1) == size(this.E_u, 1)
                m = size(u_in, 1);
                u_out = this.M_u(1:m, 1:m) * u_in;
            else
                u_out = this.M_u * u_in;
            end
        end

        function this = MD_Elliptic_u_Prior_Interface_Transient_ADR_2D(alpha_u, M, S)
            this@MD_Elliptic_u_Prior_Interface(alpha_u);
            this.M_u = M;
            this.E_u = (1.e-2) * S + M;

            num_sing_vals = 1000;
            oversampling = 20;
            num_subspace_iters = 1;
            u_vec = zeros(size(M, 1) / 2, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);
            sing_vecs_output = zeros(size(M, 1), 2 * num_sing_vals);
            sing_vals = zeros(2 * num_sing_vals, 1);
            count = 1;
            for i = 1:num_sing_vals
                sing_vecs_output(:, count) = [this.sing_vecs_output(:, i); 0 * this.sing_vecs_output(:, i)];
                sing_vals(count) = this.sing_vals(i);
                count = count + 1;
                sing_vecs_output(:, count) = [0 * this.sing_vecs_output(:, i); this.sing_vecs_output(:, i)];
                sing_vals(count) = this.sing_vals(i);
                count = count + 1;
            end
            this.sing_vecs_output = sing_vecs_output;
            this.sing_vals = sing_vals;
        end

    end

end
