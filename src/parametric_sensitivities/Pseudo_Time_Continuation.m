classdef Pseudo_Time_Continuation < handle

    properties
        psen_op
        z_nom
        theta_nom
        n
        z_k
        B_k
        grad_k
        f_k
        N
    end

    methods

        function this = Pseudo_Time_Continuation(obj, pcon, z_nom, theta_nom)
            this.psen_op = Parameteric_Sensitivity_Operators(obj, pcon);
            this.z_nom = z_nom;
            this.theta_nom = theta_nom;
            this.n = length(z_nom);
        end

        function [z_out] = Apply_Inv_Hessian(this, z_in, z, theta)

            if norm(this.psen_op.current_z - z) ~= 0 || norm(this.psen_op.current_theta - theta) ~= 0
                this.psen_op.Solve_Forward_and_Adjoint_Problems(z, theta);
            end

            tol = 1.e-7;
            max_iter = length(z);
            [z_out, flag, relres, iter, resvec] = pcg(@(x)this.psen_op.Apply_RS_Hessian(x), z_in, tol, max_iter);
            if flag ~= 0
                disp('CG did not converge');
            end
        end

        function [z] = Forward_Euler_Continuation(this, theta, N)
            this.N = N;
            this.z_k = zeros(this.n, N + 1);
            this.z_k(:, 1) = this.z_nom;

            t = linspace(0, 1, N + 1);
            theta_k = (theta - this.theta_nom) * t + this.theta_nom;
            d_theta = theta - this.theta_nom;
            dt = 1 / N;

            this.B_k = zeros(this.n, N);
            this.grad_k = zeros(length(this.z_nom), N);
            this.f_k = zeros(this.n, N);

            for k = 1:N
                this.grad_k(:, k) = this.psen_op.Solve_Forward_and_Adjoint_Problems(this.z_k(:, k), theta_k(:, k));
                this.B_k(:, k) = this.psen_op.Apply_B(d_theta);
                this.f_k(:, k) = -this.Apply_Inv_Hessian(this.B_k(:, k), this.z_k(:, k), theta_k(:, k));

                this.z_k(:, k + 1) = this.z_k(:, k) + dt * this.f_k(:, k);
            end

            z = this.z_k(:, end);
        end

    end
end
