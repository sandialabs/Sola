classdef Quasi_Newton_Preconditioner < handle

    properties
        param_qn_data
        block_qn_data
        param_current_data_step
        block_current_data_step
        tau
        max_size
        P
        W
        PW_mags
    end

    methods

        % Overload this function if a better initialization is available
        function [z_out] = Apply_Initial_Inverse_Hessian_Approximation(this, z_in)
            z_out = z_in;
        end

        function this = Quasi_Newton_Preconditioner()
            this.tau = 1.e-4;
            this.max_size = 2;
        end

        function [] = Set_N(this, N)
            this.param_qn_data = cell(N, 1);
            this.block_qn_data = cell(N, 1);
            this.param_current_data_step = 0;
            this.block_current_data_step = 0;
            this.P = [];
            this.W = [];
            this.PW_mags = [];
        end

        function [] = Add_Parametric_Quasi_Newton_Data(this, s_k, y_k)
            this.param_current_data_step = this.param_current_data_step + 1;
            this.param_qn_data{this.param_current_data_step} = struct;
            this.param_qn_data{this.param_current_data_step}.rho = 1 / (s_k' * y_k);
            if this.param_qn_data{this.param_current_data_step}.rho < 0.0
                error('Error: Negative Curvature Parameter');
            end
            this.param_qn_data{this.param_current_data_step}.s = s_k;
            this.param_qn_data{this.param_current_data_step}.y = y_k;
        end

        function [] = Add_Block_Quasi_Newton_Step(this, p, w, p_dot_w)
            num_vecs = length(this.PW_mags) + 1;
            normalization = p' * p;

            if num_vecs == 1
                this.P = p / normalization;
                this.W = w / normalization;
                this.PW_mags = p_dot_w / normalization^2;
            else
                D = zeros(num_vecs, num_vecs);
                D(1:(num_vecs - 1), 1:(num_vecs - 1)) = diag(this.PW_mags);
                D(num_vecs, 1:(num_vecs - 1)) = w' * this.P / normalization;
                D(1:(num_vecs - 1), num_vecs) = this.W' * p / normalization;
                D(num_vecs, num_vecs) = p_dot_w / normalization^2;
                [V, Lambda] = eig(D, 'vector');
                [~, J] = sort(Lambda, 'descend');
                Lambda = Lambda(J);
                V = V(:, J);
                I = find(Lambda > this.tau);
                if length(I) > this.max_size
                    I = I(1:this.max_size);
                end
                this.P = [this.P, p / normalization] * V(:, 1:I(end));
                this.W = [this.W, w / normalization] * V(:, 1:I(end));
                this.PW_mags = Lambda(1:I(end));
            end

        end

        function [] = Add_Block_Quasi_Newton_Data(this)
            this.block_current_data_step = this.block_current_data_step + 1;
            this.block_qn_data{this.block_current_data_step} = struct;
            this.block_qn_data{this.block_current_data_step}.Dr = diag(this.PW_mags);
            this.block_qn_data{this.block_current_data_step}.Pr = this.P;
            this.block_qn_data{this.block_current_data_step}.Wr = this.W;
            this.P = [];
            this.W = [];
            this.PW_mags = [];
        end

        function [z_out] = Apply_Inverse_Hessian_Approximation(this, z_in)
            [z_out] = this.Apply_QN_Inverse_Hessian_Approximation(z_in, this.param_current_data_step, this.block_current_data_step);
        end

        function [z_out] = Apply_QN_Inverse_Hessian_Approximation(this, z_in, param_counter, block_counter)
            if param_counter == 0 && block_counter == 0

                z_out = this.Apply_Initial_Inverse_Hessian_Approximation(z_in);

            elseif block_counter == param_counter

                s = this.param_qn_data{param_counter}.s;
                y = this.param_qn_data{param_counter}.y;
                rho = this.param_qn_data{param_counter}.rho;
                alpha = s(:, end)' * z_in;
                tmp = z_in - rho * alpha * y;

                tmp_out = this.Apply_QN_Inverse_Hessian_Approximation(tmp, param_counter - 1, block_counter);

                z_out = tmp_out - rho * (tmp_out' * y) * s + rho * alpha * s;

            else

                Dk = this.block_qn_data{block_counter}.Dr;
                Pk = this.block_qn_data{block_counter}.Pr;
                Wk = this.block_qn_data{block_counter}.Wr;

                tmp1 = Pk' * z_in;
                tmp1 = linsolve(Dk, tmp1);
                tmp2 = Pk * tmp1;
                tmp1 = z_in - Wk * tmp1;

                tmp_out = this.Apply_QN_Inverse_Hessian_Approximation(tmp1, param_counter, block_counter - 1);

                tmp3 = Wk' * tmp_out;
                tmp3 = linsolve(Dk, tmp3);
                tmp3 = tmp_out - Pk * tmp3;

                z_out = tmp3 + tmp2;

            end
        end

    end
end
