%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Pseudo_Time_Continuation < handle

    properties
        z_bar
        sen_op
        qn_prec
        use_qn_prec
        print_output
        print_iter
    end

    methods

        function this = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec, print_output, print_iter)
            this.z_bar = z_bar;
            this.sen_op = sen_op;
            this.qn_prec = qn_prec;
            this.use_qn_prec = true;
            if nargin < 4
                this.print_output = false;
            else
                this.print_output = print_output;
            end
            if nargin < 5
                this.print_iter = false;
            else
                this.print_iter = print_iter;
            end
        end

        function [z_out] = Apply_Inverse_Hessian(this, z_in, z, theta_traj, time_index)

            tol = 1.e-6;
            max_iter = length(z_in);

            z_out = 0 * z_in;
            r = z_in;
            v = this.qn_prec.Apply_Inverse_Hessian_Approximation(r);
            p = v;
            scalar = r' * v;
            rel_tol = tol * norm(z_in);
            iter = 0;

            while (sqrt(scalar) > rel_tol) && (norm(r) > rel_tol) && (iter < max_iter)
                iter = iter + 1;
                if this.print_iter
                    disp(['Iteration = ', num2str(iter)]);
                    disp(['Sqrt(rho) = ', num2str(sqrt(scalar))]);
                    disp(['Norm(r) = ', num2str(norm(r))]);
                end

                w = this.sen_op.Apply_Hessian(p, z, theta_traj, time_index);

                p_dot_w = w' * p;
                this.qn_prec.Add_Block_Quasi_Newton_Step(p, w, p_dot_w);
                alpha = scalar / p_dot_w;
                z_out = z_out + alpha * p;
                r = r - alpha * w;
                v = this.qn_prec.Apply_Inverse_Hessian_Approximation(r);
                scalar_old = scalar;
                scalar = v' * r;
                p = v + (scalar / scalar_old) * p;
            end

            if this.print_output
                disp(['Total iterations = ', num2str(iter)]);
                disp(['Relative Residual = ', num2str(norm(r) / norm(z_in))]);
            end
        end

        function [z_k, grad_k] = Pseudo_Time_Continuation_Forward_Euler(this, theta_traj)

            N = theta_traj.Get_Number_of_Timesteps();
            this.qn_prec.Set_N(N);

            z_k = zeros(length(this.z_bar), N + 1);
            z_k(:, 1) = this.z_bar;

            grad_k = zeros(length(this.z_bar), N + 1);
            dt = 1 / N;

            grad_k(:, 1) = this.sen_op.Gradient(z_k(:, 1), theta_traj, 0);
            for k = 1:N

                B = this.sen_op.Apply_B(z_k(:, k), theta_traj, k - 1);
                f = this.Apply_Inverse_Hessian(B, z_k(:, k), theta_traj, k - 1);

                z_tmp = z_k(:, k) - dt * f;
                grad_tmp = this.sen_op.Gradient(z_tmp, theta_traj, k);

                if this.use_qn_prec
                    this.qn_prec.Add_Block_Quasi_Newton_Data();
                    s_k = z_tmp - z_k(:, k);
                    y_k = grad_tmp - grad_k(:, k) - dt * B;
                    this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                end

                f = this.Apply_Inverse_Hessian(grad_tmp, z_tmp, theta_traj, k);
                z_k(:, k + 1) = z_tmp - f;
                grad_k(:, k + 1) = this.sen_op.Gradient(z_k(:, k + 1), theta_traj, k);

                if this.use_qn_prec
                    this.qn_prec.Add_Block_Quasi_Newton_Data();
                    s_k = z_k(:, k + 1) - z_tmp;
                    y_k = grad_k(:, k + 1) - grad_tmp;
                    this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                end

            end
        end

        function [z_k, grad_k] = Pseudo_Time_Continuation_Modified_Euler(this, theta_traj)

            N = theta_traj.Get_Number_of_Timesteps();
            this.qn_prec.Set_N(N);

            z_k = zeros(length(this.z_bar), N + 1);
            z_k(:, 1) = this.z_bar;

            grad_k = zeros(length(this.z_bar), N + 1);
            dt = 1 / N;

            grad_k(:, 1) = this.sen_op.Gradient(z_k(:, 1), theta_traj, 0);
            for k = 1:N

                B = this.sen_op.Apply_B(z_k(:, k), theta_traj, k - 1);
                f = this.Apply_Inverse_Hessian(B, z_k(:, k), theta_traj, k - 1);
                z_tmp1 = z_k(:, k) - 0.5 * dt * f;

                grad_tmp1 = this.sen_op.Gradient(z_tmp1, theta_traj, k - .5);

                if this.use_qn_prec
                    this.qn_prec.Add_Block_Quasi_Newton_Data();
                    s_k = z_tmp1 - z_k(:, k);
                    y_k = grad_tmp1 - grad_k(:, k) - 0.5 * dt * B;
                    this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                end

                B = this.sen_op.Apply_B(z_tmp1, theta_traj, k - .5);
                f = this.Apply_Inverse_Hessian(B, z_tmp1, theta_traj, k - .5);
                z_tmp2 = z_k(:, k) - dt * f;

                grad_tmp2 = this.sen_op.Gradient(z_tmp2, theta_traj, k);

                if this.use_qn_prec
                    this.qn_prec.Add_Block_Quasi_Newton_Data();
                    s_k = z_tmp2 - z_tmp1;
                    y_k = grad_tmp2 - grad_tmp1 - 0.5 * dt * B;
                    this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                end

                f = this.Apply_Inverse_Hessian(grad_tmp2, z_tmp2, theta_traj, k);
                z_k(:, k + 1) = z_tmp2 - f;
                grad_k(:, k + 1) = this.sen_op.Gradient(z_k(:, k + 1), theta_traj, k);

                if this.use_qn_prec
                    this.qn_prec.Add_Block_Quasi_Newton_Data();
                    s_k = z_k(:, k + 1) - z_tmp2;
                    y_k = grad_k(:, k + 1) - grad_tmp2;
                    this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                end

            end
        end

    end
end
