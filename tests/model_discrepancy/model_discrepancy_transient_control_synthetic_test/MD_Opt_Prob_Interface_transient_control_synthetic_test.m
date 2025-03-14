classdef MD_Opt_Prob_Interface_transient_control_synthetic_test < MD_Opt_Prob_Interface

    properties
        n_y
        n_t
        x
        t
        J
        T
    end

    methods (Access = public)

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            z_out = this.J' * u_in;
        end

        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            z_out = this.J' * this.J * z_in;
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            grad_u = u - this.T;
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = u_in;
        end

    end

    methods

        function this = MD_Opt_Prob_Interface_transient_control_synthetic_test(n_y,n_t)
            this.n_y = n_y;
            this.x = linspace(0, 1, n_y)';
            this.n_t = n_t;
            this.t = linspace(0,1,n_t)';

            J1 = 1 - (0:(n_y-1))/(n_y-1);
            J2 = (0:(n_y-1))/(n_y-1);
            this.J = kron(eye(n_t),[J1;J2])';

            tmp = 2*J2'*this.t' + ones(n_y,1)*this.t';
            this.T = tmp(:);
        end

    end

end
