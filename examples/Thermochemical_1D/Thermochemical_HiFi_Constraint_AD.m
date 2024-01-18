classdef Thermochemical_HiFi_Constraint_AD < Dynamic_Constraint_AD

    properties
        fe
        control_time_nodes

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

        legendre_polys
    end

    methods (Access = public)

        function [f] = Time_Instance_RHS_AD(this, y, z, t)
            T = this.I_T * y;
            u1 = this.I_u1 * y;
            u2 = this.I_u2 * y;
            v1 = this.I_v1 * y;
            v2 = this.I_v2 * y;
            v3 = this.I_v3 * y;

            w = this.Temporal_Weights(t);
            zt = reshape(z, this.fe.m, this.control_time_nodes) * w;
            heating = this.fe.M * zt;

            T_coll = this.fe.nodes_to_coll_points * T;
            v1_coll = this.fe.nodes_to_coll_points * v1;
            v2_coll = this.fe.nodes_to_coll_points * v2;
            v3_coll = this.fe.nodes_to_coll_points * v3;

            R_coll = this.Evaluate_Arrhenius_Law(T_coll);

            reaction_1 = this.react_rate_1 * this.fe.coll_point_integration * (R_coll .* v1_coll .* v2_coll);
            reaction_2 = this.react_rate_2 * this.fe.coll_point_integration * (R_coll .* v1_coll .* v3_coll);

            rhs_T = -this.diff_T * this.fe.S * T + heating - this.cooling_1 * reaction_1 - this.cooling_2 * reaction_2;
            rhs_u1 = -this.diff_u1 * this.fe.S * u1 + reaction_1;
            rhs_u2 = -this.diff_u2 * this.fe.S * u2 + reaction_2;
            rhs_v1 = -this.diff_v1 * this.fe.S * v1 - reaction_1 - reaction_2;
            rhs_v2 = -this.diff_v2 * this.fe.S * v2 - reaction_1;
            rhs_v3 = -this.diff_v3 * this.fe.S * v3 - reaction_2;

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
            v1_0 = 10 * ones(this.fe.m, 1) + 30 * exp(-80 * (this.fe.x - 0.3).^2) + 30 * exp(-70 * (this.fe.x - 0.6).^2);
        end

        function [v2_0] = Initial_v2(this)
            v2_0 = 10 * ones(this.fe.m, 1) + 40 * exp(-80 * (this.fe.x - 0.3).^2) + 40 * exp(-70 * (this.fe.x - 0.6).^2);
        end

        function [v3_0] = Initial_v3(this)
            v3_0 = 10 * ones(this.fe.m, 1) + 2 * exp(-80 * (this.fe.x - 0.3).^2) + 90 * exp(-70 * (this.fe.x - 0.6).^2);
        end

    end

    methods (Access = public)

        function this = Thermochemical_HiFi_Constraint_AD(n_y, control_time_nodes, T, n_t)
            this@Dynamic_Constraint_AD(6 * n_y, n_y * control_time_nodes, T, n_t);
            this.control_time_nodes = control_time_nodes;

            this.fe = Linear_1D_Finite_Elements(n_y);

            this.I_T = zeros(n_y, 6 * n_y);
            this.I_T(:, 1:n_y) = eye(n_y);

            this.I_u1 = zeros(n_y, 6 * n_y);
            this.I_u1(:, (n_y + 1):(2 * n_y)) = eye(n_y);

            this.I_u2 = zeros(n_y, 6 * n_y);
            this.I_u2(:, (2 * n_y + 1):(3 * n_y)) = eye(n_y);

            this.I_v1 = zeros(n_y, 6 * n_y);
            this.I_v1(:, (3 * n_y + 1):(4 * n_y)) = eye(n_y);

            this.I_v2 = zeros(n_y, 6 * n_y);
            this.I_v2(:, (4 * n_y + 1):(5 * n_y)) = eye(n_y);

            this.I_v3 = zeros(n_y, 6 * n_y);
            this.I_v3(:, (5 * n_y + 1):(6 * n_y)) = eye(n_y);

            this.arrhenius_scale = 10.0;

            this.diff_T = 1.e-4;
            this.diff_u1 = 1.e-5;
            this.diff_u2 = 1.e-5;
            this.diff_v1 = 1.e-5;
            this.diff_v2 = 1.e-5;
            this.diff_v3 = 1.e-5;

            this.cooling_1 = 0.5;
            this.cooling_2 = 1.0;

            this.react_rate_1 = 5.e-2;
            this.react_rate_2 = 1.e-1;

            legendre_polys = zeros(this.n_t, this.control_time_nodes);
            for k = 1:this.control_time_nodes
                legendre_polys(:, k) = legendreP(k - 1, 2 * this.t_mesh - 1);
            end
            this.legendre_polys = legendre_polys;

        end

        function [f] = Map_Controller_to_Mesh(this, z)
            zt = reshape(z, this.fe.m, this.control_time_nodes);
            f = zt * this.legendre_polys';
        end

        function [w] = Temporal_Weights(this, t)
            w = this.Legendre_Polynomials(t);
        end

        function [P] = Legendre_Polynomials(this, t)
            t = 2 * t - 1;
            P = [0 * t + 1.0; t];
            for k = 3:this.control_time_nodes
                n = k - 1;
                P = [P; (2 * n - 1) / n * t * P(k - 1) - (n - 1) / n * P(k - 2)];
            end
        end

    end

end
