classdef Pseudo_Time_Continuation_BFGS < Pseudo_Time_Continuation

    properties
        rho_k
        s_k
        y_k
        iter_k
        recursion_k
    end

    methods (Abstract, Access = public)

        [z_out] = Apply_Nominal_Inv_Hessian(this, z_in)

    end

    methods

        function this = Pseudo_Time_Continuation_BFGS(obj, pcon, z_nom, theta_nom)
            this@Pseudo_Time_Continuation(obj, pcon, z_nom, theta_nom);
            this.iter_k = 1;
        end

        function [z_out] = Apply_Inv_Hessian(this, z_in, z, theta)

            if this.iter_k == 1
                this.rho_k = zeros(this.N, 1);
                this.s_k = zeros(this.n, this.N);
                this.y_k = zeros(this.n, this.N);
                this.recursion_k = 1;
            end

            if this.recursion_k == 1
                z_out = this.Apply_Nominal_Inv_Hessian(z_in);
                this.iter_k = this.iter_k + 1;
                this.recursion_k = this.iter_k;
            else
                k = this.recursion_k;
                dt = 1 / this.N;
                this.s_k(:, k - 1) = this.z_k(:, k) - this.z_k(:, k - 1);
                this.y_k(:, k - 1) = this.grad_k(:, k) - this.grad_k(:, k - 1) - dt * this.B_k(:, k - 1);
                this.rho_k(k - 1) = 1 / (this.s_k(:, k - 1)' * this.y_k(:, k - 1));
                if ~(this.rho_k(k - 1) > 0)
                    error(['Error: rho = ', num2str(this.rho_k(k - 1))]);
                end

                alpha = this.s_k(:, k - 1)' * z_in;
                tmp = z_in - this.rho_k(k - 1) * alpha * this.y_k(:, k - 1);
                this.recursion_k = this.recursion_k - 1;
                tmp_out = this.Apply_Inv_Hessian(tmp);
                z_out = tmp_out - this.rho_k(k - 1) * (tmp_out' * this.y_k(:, k - 1)) * this.s_k(:, k - 1) + this.rho_k(k - 1) * alpha * this.s_k(:, k - 1);
            end
        end

    end
end
