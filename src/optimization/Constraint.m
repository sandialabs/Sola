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
        [u] = State_Solve(this, z);

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % Mv: c_u(u, z)^{-T}v in R^{n_u}
        [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z);

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % Mv: c_z(u, z)^{T}v in R^{n_z}
        [Mv] = c_z_Transpose_Apply(this, v, u, z);

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % Mv: c_u(u, z)^{-1}v in R^{n_u}
        [Mv] = c_u_Inverse_Apply(this, v, u, z);

        % Input:
        % v: a direction v in R^{n_z}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % Output:
        % Mv: c_z(u, z)v in R^{n_u}
        [Mv] = c_z_Apply(this, v, u, z);

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % lambda: the adjoint state lambda in R^{n_u}
        % Output:
        % Mv: lambda^T c_{u, u}(u, z)v in R^{n_u}
        [Mv] = c_uu_Apply(this, v, u, z, lambda);

        % Input:
        % v: a direction v in R^{n_z}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % lambda: the adjoint state lambda in R^{n_u}
        % Output:
        % Mv: lambda^T c_{u, z}(u, z)v in R^{n_u}
        [Mv] = c_uz_Apply(this, v, u, z, lambda);

        % Input:
        % v: a direction v in R^{n_u}
        % u: the state u in R^{n_u}
        % z: the control z in R^{n_z}
        % lambda: the adjoint state lambda in R^{n_u}
        % Output:
        % Mv: lambda^T c_{z, u}(u, z)v in R^{n_z}
        [Mv] = c_zu_Apply(this, v, u, z, lambda);

        % Input:
        % * v: a direction v in R^{n_z}
        % * u: the state u in R^{n_u}
        % * z: the control z in R^{n_z}
        % * lambda: the adjoint state lambda in R^{n_u}
        % Output:
        % Mv: lambda^T c_{z, z}(u, z)v in R^{n_z}
        [Mv] = c_zz_Apply(this, v, u, z, lambda);

    end

    methods (Access = public)

        function this = Constraint( )

        end

    end
end
