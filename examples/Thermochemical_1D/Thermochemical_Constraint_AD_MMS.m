classdef Thermochemical_Constraint_AD_MMS < Dynamic_Constraint_AD

    properties
        fe

        I_T
        I_u1
        I_u2
        I_v1
        I_v2
        I_v3

        arrhenius_scale

        diff_T
        diff_u1
        diff_u2
        diff_v1
        diff_v2
        diff_v3

        cooling_1
        cooling_2

        react_rate_1
        react_rate_2
    end

    methods (Access = public)

        function [f] = Time_Instance_RHS_AD(this, y, z, t)
            T = this.I_T * y;
            u1 = this.I_u1 * y;
            u2 = this.I_u2 * y;
            v1 = this.I_v1 * y;
            v2 = this.I_v2 * y;
            v3 = this.I_v3 * y;

            z_T = kron(eye(this.n_t), this.I_T) * z;
            z_u1 = kron(eye(this.n_t), this.I_u1) * z;
            z_u2 = kron(eye(this.n_t), this.I_u2) * z;
            z_v1 = kron(eye(this.n_t), this.I_v1) * z;
            z_v2 = kron(eye(this.n_t), this.I_v2) * z;
            z_v3 = kron(eye(this.n_t), this.I_v3) * z;

            w = this.Temporal_Weights(t);

            zt_T = reshape(z_T, this.fe.m, this.n_t) * w;
            heating = this.fe.M * zt_T;

            zt_u1 = reshape(z_u1, this.fe.m, this.n_t) * w;
            u1_forcing = this.fe.M * zt_u1;

            zt_u2 = reshape(z_u2, this.fe.m, this.n_t) * w;
            u2_forcing = this.fe.M * zt_u2;

            zt_v1 = reshape(z_v1, this.fe.m, this.n_t) * w;
            v1_forcing = this.fe.M * zt_v1;

            zt_v2 = reshape(z_v2, this.fe.m, this.n_t) * w;
            v2_forcing = this.fe.M * zt_v2;

            zt_v3 = reshape(z_v3, this.fe.m, this.n_t) * w;
            v3_forcing = this.fe.M * zt_v3;

            T_coll = this.fe.nodes_to_coll_points * T;
            v1_coll = this.fe.nodes_to_coll_points * v1;
            v2_coll = this.fe.nodes_to_coll_points * v2;
            v3_coll = this.fe.nodes_to_coll_points * v3;

            R_coll = this.Evaluate_Arrhenius_Law(T_coll);

            reaction_1 = this.react_rate_1 * this.fe.coll_point_integration * (R_coll .* v1_coll .* v2_coll);
            reaction_2 = this.react_rate_2 * this.fe.coll_point_integration * (R_coll .* v1_coll .* v3_coll);

            rhs_T = -this.diff_T * this.fe.S * T + heating - this.cooling_1 * reaction_1 - this.cooling_2 * reaction_2;
            rhs_u1 = -this.diff_u1 * this.fe.S * u1 + reaction_1 + u1_forcing;
            rhs_u2 = -this.diff_u2 * this.fe.S * u2 + reaction_2 + u2_forcing;
            rhs_v1 = -this.diff_v1 * this.fe.S * v1 - reaction_1 - reaction_2 + v1_forcing;
            rhs_v2 = -this.diff_v2 * this.fe.S * v2 - reaction_1 + v2_forcing;
            rhs_v3 = -this.diff_v3 * this.fe.S * v3 - reaction_2 + v3_forcing;

            f_T = this.fe.M \ rhs_T;
            f_u1 = this.fe.M \ rhs_u1;
            f_u2 = this.fe.M \ rhs_u2;
            f_v1 = this.fe.M \ rhs_v1;
            f_v2 = this.fe.M \ rhs_v2;
            f_v3 = this.fe.M \ rhs_v3;

            f = this.I_T' * f_T + this.I_u1' * f_u1 + this.I_u2' * f_u2 + this.I_v1' * f_v1 + this.I_v2' * f_v2 + this.I_v3' * f_v3;
        end

        function [h] = Initial_Condition_AD(this, z)
            T_0 = this.Initial_T();
            u1_0 = this.Initial_u1();
            u2_0 = this.Initial_u2();
            v1_0 = this.Initial_v1();
            v2_0 = this.Initial_v2();
            v3_0 = this.Initial_v3();

            h = this.I_T' * T_0 + this.I_u1' * u1_0 + this.I_u2' * u2_0 + this.I_v1' * v1_0 + this.I_v2' * v2_0 + this.I_v3' * v3_0;
        end

        function [R] = Evaluate_Arrhenius_Law(this, T)
            R = exp(-this.arrhenius_scale ./ T);
        end

        function [T_0] = Initial_T(this)
            T_0 = ones(this.fe.m, 1);
        end

        function [u1_0] = Initial_u1(this)
            u1_0 = zeros(this.fe.m, 1);
        end

        function [u2_0] = Initial_u2(this)
            u2_0 = zeros(this.fe.m, 1);
        end

        function [v1_0] = Initial_v1(this)
            v1_0 = 20 * ones(this.fe.m, 1);
        end

        function [v2_0] = Initial_v2(this)
            v2_0 = 10 * ones(this.fe.m, 1);
        end

        function [v3_0] = Initial_v3(this)
            v3_0 = 12 * ones(this.fe.m, 1);
        end

    end

    methods (Access = public)

        function this = Thermochemical_Constraint_AD_MMS(n_y, n_z, T, n_t)
            this@Dynamic_Constraint_AD(6 * n_y, n_z, T, n_t);

            this.fe = Linear_1D_Finite_Elements(n_y);

            this.I_T = sparse(n_y, 6 * n_y);
            this.I_T(:, 1:n_y) = eye(n_y);

            this.I_u1 = sparse(n_y, 6 * n_y);
            this.I_u1(:, (n_y + 1):(2 * n_y)) = eye(n_y);

            this.I_u2 = sparse(n_y, 6 * n_y);
            this.I_u2(:, (2 * n_y + 1):(3 * n_y)) = eye(n_y);

            this.I_v1 = sparse(n_y, 6 * n_y);
            this.I_v1(:, (3 * n_y + 1):(4 * n_y)) = eye(n_y);

            this.I_v2 = sparse(n_y, 6 * n_y);
            this.I_v2(:, (4 * n_y + 1):(5 * n_y)) = eye(n_y);

            this.I_v3 = sparse(n_y, 6 * n_y);
            this.I_v3(:, (5 * n_y + 1):(6 * n_y)) = eye(n_y);

            this.arrhenius_scale = .1;

            this.diff_T = 1.13;
            this.diff_u1 = .19;
            this.diff_u2 = 3.1;
            this.diff_v1 = 4.0;
            this.diff_v2 = 8.3;
            this.diff_v3 = .5;

            this.cooling_1 = .1;
            this.cooling_2 = 2;

            this.react_rate_1 = .44;
            this.react_rate_2 = 5.4;

        end

        function [w] = Temporal_Weights(this, t)
            tmp = (this.t_mesh - t) / (this.t_mesh(2) - this.t_mesh(1));
            w = 1 - min(abs(tmp), 1);
        end

    end

end
