classdef MD_Opt_Prob_Interface_transient_multi_state_synthetic_test < MD_Opt_Prob_Interface

    properties
        n_y
        n_t
        x
        M
        c_low
    end

    methods (Access = public)

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            J = zeros(2*this.n_y*this.n_t,this.n_y);
            for k = 1:this.n_t
                I = (1:this.n_y) + 2*(k-1)*this.n_y;
                J(I,:) = (this.c_low^(k-1)) * 3 * diag(z.^2);
                I = ((this.n_y+1):(2*this.n_y)) + 2*(k-1)*this.n_y;
                J(I,:) = (this.c_low^k) * 3 * diag(z.^2);
            end
            z_out = J' * u_in;
        end

        % This implementation assumes that it is evaluated at the optimal z so that
        % the adjoint=0, a more general impl a term multiplied by the adjoint variable
        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            tmp1 = (this.c_low^this.n_t) * 3 * diag(z.^2) * z_in;
            tmp2 = this.M * tmp1;
            z_out = (this.c_low^this.n_t) * 3 * diag(z.^2) * tmp2;
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            grad_u = 0 * u;
            grad_u((end-this.n_y+1):end) = this.M * (u((end-this.n_y+1):end) - (this.c_low^this.n_t)*(1 + this.x).^3);
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = 0 * u_in;
            u_out((end-this.n_y+1):end,:) = this.M * u_in((end-this.n_y+1):end,:);
        end

    end

    methods

        function this = MD_Opt_Prob_Interface_transient_multi_state_synthetic_test(n_y,n_t,x,M,c_low)
            this.n_y = n_y;
            this.n_t = n_t;
            this.x = x;
            this.M = M;
            this.c_low = c_low;
        end

    end

end
