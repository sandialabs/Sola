classdef Pseudo_Time_Continuation < handle

    properties
        z_bar
        theta_bar
        sen_op
        use_bfgs_prec
    end

    methods

        % Overload this function if a better initialization is available
        function [z_out] = Apply_Initial_Inverse_BFGS_Hessian(this, z_in)
            z_out = z_in;
        end

        function this = Pseudo_Time_Continuation(z_bar, theta_bar, sen_op)
            this.z_bar = z_bar;
            this.theta_bar = theta_bar;
            this.sen_op = sen_op;
            this.use_bfgs_prec = true;
        end

        function [z_out] = Apply_Inverse_Hessian(this, z_in, z, theta, rho_k, s_k, y_k)

            tol = 1.e-5;
            max_iter = length(z_in);

            print_iter = false;
            print_output = false;

            z_out = 0 * z_in;
            r = z_in;
            v = this.Apply_Inverse_BFGS_Hessian(r, rho_k, s_k, y_k);
            p = v;
            scalar = r' * v;
            rel_tol = tol * norm(z_in);
            iter = 0;
            while (sqrt(scalar) > rel_tol) && (norm(r) > rel_tol) && (iter < max_iter)
                iter = iter + 1;
                if print_iter
                    disp(['Iteration = ', num2str(iter)]);
                    disp(['Sqrt(rho) = ', num2str(sqrt(scalar))]);
                    disp(['Norm(r) = ', num2str(norm(r))]);
                end

                w = this.sen_op.Apply_Hessian(p, z, theta);
                alpha = scalar / (w' * p);
                z_out = z_out + alpha * p;
                r = r - alpha * w;
                v = this.Apply_Inverse_BFGS_Hessian(r, rho_k, s_k, y_k);
                scalar_old = scalar;
                scalar = v' * r;
                p = v + (scalar / scalar_old) * p;
            end
            if print_output
                disp(['Total iterations = ', num2str(iter)]);
                disp(['Relative Residual = ', num2str(sqrt(scalar) / norm(z_in))]);
            end
        end

        function [z_k, grad_k] = Pseudo_Time_Continuation_Forward_Euler(this, theta_star, N)
            z_k = zeros(length(this.z_bar), N + 1);
            z_k(:, 1) = this.z_bar;

            grad_k = zeros(length(this.z_bar), N + 1);
            rho_k = zeros(N, 1);
            s_k = zeros(length(this.z_bar), N);
            y_k = zeros(length(this.z_bar), N);

            d_theta = theta_star - this.theta_bar;
            dt = 1 / N;

            grad_k(:, 1) = this.sen_op.Gradient(z_k(:, 1), this.theta_bar);
            for k = 1:N

                B = this.sen_op.Apply_B(d_theta, z_k(:, k), this.theta_bar + (k - 1) * dt * d_theta);
                f = -this.Apply_Inverse_Hessian(B, z_k(:, k), this.theta_bar + (k - 1) * dt * d_theta, rho_k(1:(k - 1)), s_k(:, 1:(k - 1)), y_k(:, 1:(k - 1)));

                z_k(:, k + 1) = z_k(:, k) + dt * f;
                grad_k(:, k + 1) = this.sen_op.Gradient(z_k(:, k + 1), this.theta_bar + k * dt * d_theta);

                s_k(:, k) = z_k(:, k + 1) - z_k(:, k);
                y_k(:, k) = grad_k(:, k + 1) - grad_k(:, k) - dt * B;
                rho_k(k) = 1 / (s_k(:, k)' * y_k(:, k));
                if rho_k(k) < 0.0
                    error(['Error: rho = ', num2str(rho_k(k))]);
                end
            end
        end

        function [z_k, grad_k] = Pseudo_Time_Continuation_Modified_Euler(this, theta_star, N)
            z_k = zeros(length(this.z_bar), N + 1);
            z_k(:, 1) = this.z_bar;

            grad_k = zeros(length(this.z_bar), N + 1);
            rho_k = zeros(N, 1);
            s_k = zeros(length(this.z_bar), N);
            y_k = zeros(length(this.z_bar), N);

            d_theta = theta_star - this.theta_bar;
            dt = 1 / N;

            grad_k(:, 1) = this.sen_op.Gradient(z_k(:, 1), this.theta_bar);
            for k = 1:N

                B_tmp = this.sen_op.Apply_B(d_theta, z_k(:, k), this.theta_bar + (k - 1) * dt * d_theta);
                f_tmp = -this.Apply_Inverse_Hessian(B_tmp, z_k(:, k), this.theta_bar + (k - 1) * dt * d_theta, rho_k(1:(k - 1)), s_k(:, 1:(k - 1)), y_k(:, 1:(k - 1)));
                z_tmp = z_k(:, k) + 0.5 * dt * f_tmp;

                this.sen_op.Gradient(z_tmp, this.theta_bar + (k - .5) * dt * d_theta);
                B = this.sen_op.Apply_B(d_theta, z_tmp, this.theta_bar + (k - .5) * dt * d_theta);
                f = -this.Apply_Inverse_Hessian(B, z_tmp, this.theta_bar + (k - .5) * dt * d_theta, rho_k(1:(k - 1)), s_k(:, 1:(k - 1)), y_k(:, 1:(k - 1)));
                z_k(:, k + 1) = z_k(:, k) + dt * f;

                grad_k(:, k + 1) = this.sen_op.Gradient(z_k(:, k + 1), this.theta_bar + k * dt * d_theta);

                s_k(:, k) = z_k(:, k + 1) - z_k(:, k);
                % I'm seeing better imperical results using "B" rather than
                % "B_tmp" but the theory indicates that "B_tmp" is the
                % correct term. I need to look into this more deeply.
                y_k(:, k) = grad_k(:, k + 1) - grad_k(:, k) - dt * B_tmp;
                rho_k(k) = 1 / (s_k(:, k)' * y_k(:, k));
                if rho_k(k) < 0.0
                    error(['Error: rho = ', num2str(rho_k(k))]);
                end
            end
        end

        function [z_out] = Apply_Inverse_BFGS_Hessian(this, z_in, rho, s, y)
            if isempty(rho) || (~this.use_bfgs_prec)
                z_out = this.Apply_Initial_Inverse_BFGS_Hessian(z_in);
            else
                alpha = s(:, end)' * z_in;
                tmp = z_in - rho(end) * alpha * y(:, end);
                if length(rho) > 1
                    tmp_out = this.Apply_Inverse_BFGS_Hessian(tmp, rho(1:(end - 1)), s(:, 1:(end - 1)), y(:, 1:(end - 1)));
                else
                    tmp_out = this.Apply_Initial_Inverse_BFGS_Hessian(tmp);
                end
                z_out = tmp_out - rho(end) * (tmp_out' * y(:, end)) * s(:, end) + rho(end) * alpha * s(:, end);
            end
        end

    end
end
