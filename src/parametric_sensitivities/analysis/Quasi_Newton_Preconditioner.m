classdef Quasi_Newton_Preconditioner < handle

    properties
        param_qn_data
        block_qn_data
        param_current_data_step
        block_current_data_step
        tau
    end

    methods

        % Overload this function if a better initialization is available
        function [z_out] = Apply_Initial_Inverse_Hessian_Approximation(this, z_in)
            z_out = z_in;
        end

        function this = Quasi_Newton_Preconditioner()
            this.tau = 1.e-6;
        end

        function [] = Set_N(this, N)
            this.param_qn_data = cell(N, 1);
            this.block_qn_data = cell(N, 1);
            this.param_current_data_step = 0;
            this.block_current_data_step = 0;
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

        function [] = Add_Block_Quasi_Newton_Data(this, P, W)

            A = P' * W;
            [L, D, perm] = ldl(A, 'vector');
            tmp = vecnorm(P)';
            tmp = tmp(perm);
            indices = diag(D) > this.tau * tmp;
            Pr = P(:, perm(indices));
            Wr = W(:, perm(indices));
            Lr = L(indices, indices);
            Dr = D(indices, indices);

            this.block_current_data_step = this.block_current_data_step + 1;
            this.block_qn_data{this.block_current_data_step} = struct;
            this.block_qn_data{this.block_current_data_step}.Lr = Lr;
            this.block_qn_data{this.block_current_data_step}.Dr = Dr;
            this.block_qn_data{this.block_current_data_step}.Pr = Pr;
            this.block_qn_data{this.block_current_data_step}.Wr = Wr;
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

                L = this.block_qn_data{block_counter}.Lr;
                D = this.block_qn_data{block_counter}.Dr;
                P = this.block_qn_data{block_counter}.Pr;
                W = this.block_qn_data{block_counter}.Wr;

                tmp1 = linsolve(L, P' * z_in);
                tmp1 = linsolve(D, tmp1);
                tmp1 = linsolve(L', tmp1);
                tmp2 = P * tmp1;
                tmp1 = z_in - W * tmp1;

                tmp_out = this.Apply_QN_Inverse_Hessian_Approximation(tmp1, param_counter, block_counter - 1);

                tmp3 = linsolve(L, W' * tmp_out);
                tmp3 = linsolve(D, tmp3);
                tmp3 = linsolve(L', tmp3);
                tmp3 = tmp_out - P * tmp3;

                z_out = tmp3 + tmp2;

            end
        end

    end
end
