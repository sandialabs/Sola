% Define the constraint function c(u, z) where
% u in R^{n_u}
% z in R^{n_z}
% c(u, z) in R^{n_u}

classdef Constraint < handle

    methods (Abstract, Access = public)

        % Input:
        % z: the control z in R^{n_z}
        % Output:
        % u: u = S(z) in R^{n_u}
        [u] = State_Solve(this, z)

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % Mv: c_u(u, z)^{-T}v in R^{n_u}
        [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % Mv: c_z(u, z)^{T}v in R^{n_z}
        [Mv] = c_z_Transpose_Apply(this, v, u, z)

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % Mv: c_u(u, z)^{-1}v in R^{n_u}
        [Mv] = c_u_Inverse_Apply(this, v, u, z)

        % Input:
        % v: a direction v in R^{n_z}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % Mv: c_z(u, z)v in R^{n_u}
        [Mv] = c_z_Apply(this, v, u, z)

    end

    methods (Access = public)

        % Input:
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % c: c(u, z)v in R^{n_u}
        function [c] = c(this, u, z)
            disp('Error: c is not implemented');
            c = 'Not Implemented';
        end

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % lambda: the adjoint state lambda in R^{n_u}
        % Output:
        % Mv: lambda^T c_{u, u}(u, z)v in R^{n_u}
        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            disp('Error: c_uu is not implemented');
            Mv = 'Not Implemented';
        end

        % Input:
        % v: a direction v in R^{n_z}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % lambda: the adjoint state lambda in R^{n_u}
        % Output:
        % Mv: lambda^T c_{u, z}(u, z)v in R^{n_u}
        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            disp('Error: c_uz is not implemented');
            Mv = 'Not Implemented';
        end

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % lambda: the adjoint state lambda in R^{n_u}
        % Output:
        % Mv: lambda^T c_{z, u}(u, z)v in R^{n_z}
        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            disp('Error: c_zu is not implemented');
            Mv = 'Not Implemented';
        end

        % Input:
        % * v: a direction v in R^{n_z}
        % * u: the state u in R^{n_u}
        % * z: the control z in R^{n_z}
        % * lambda: the adjoint state lambda in R^{n_u}
        % Output:
        % Mv: lambda^T c_{z, z}(u, z)v in R^{n_z}
        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            disp('Error: c_zz is not implemented');
            Mv = 'Not Implemented';
        end

    end

    methods (Access = public)

        function this = Constraint()

        end

        function [diffs_z, jacobian_z_transpose_check, diffs_u, jacobian_u_transpose_check, solve_res] = Finite_Difference_Constraint_Check(this, u, z)

            % c_z_Apply check
            c = this.c(u, z);
            h = -2:-1:-6;
            p = length(h);
            v = randn(length(z), 1);
            cz_v = this.c_z_Apply(v, u, z);
            diffs_z = zeros(p, 1);
            for k = 1:p
                ck =  this.c(u, z + h(k) * v);
                diffs_z(k) = norm(cz_v - (ck - c) / h(k)) / norm(cz_v);
            end
            disp('Constraint z Jacobian finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_z(k))]);
            end
            disp(' ');

            % c_z_Transpose_Apply check
            vu = randn(length(u), 1);
            vu_cz_v = vu' * cz_v;
            ctz_vu = this.c_z_Transpose_Apply(vu, u, z);
            v_ctz_vu = v' * ctz_vu;
            jacobian_z_transpose_check = abs(vu_cz_v - v_ctz_vu) / abs(vu_cz_v);
            if jacobian_z_transpose_check > 1.e-13
                disp('Error in c_z_Transpose_Apply');
            end

            % c_u_Inverse_Apply check
            diffs_u = zeros(p, 1);
            for k = 1:p
                cu_vu_fd = (this.c(u + h(k) * vu, z) - c) / h(k);
                diffs_u(k) = norm(this.c_u_Inverse_Apply(cu_vu_fd, u, z) - vu) / norm(vu);
            end
            disp('Constraint u Jacobian Inverse finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_u(k))]);
            end
            disp(' ');

            % c_u_Transpose_Inverse_Apply check
            cu_inv_vu = this.c_u_Inverse_Apply(vu, u, z);
            d = randn(length(u), 1);
            cu_trans_inv_d = this.c_u_Transpose_Inverse_Apply(d, u, z);
            jacobian_u_transpose_check = abs(vu' * cu_trans_inv_d - d' * cu_inv_vu) / abs(d' * cu_inv_vu);
            if jacobian_u_transpose_check > 1.e-13
                disp('Error in c_u_Transpose_Inverse_Apply');
            end

            % State solve check
            u_solve = this.State_Solve(z);
            solve_res = norm(this.c(u_solve, z));
            if solve_res > 1.e-10
                disp('Error in state solve');
            end
        end

    end
end
