%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Thermochemical_LoFi_Constraint_AD < Dynamic_Constraint_AD

    properties
        con_hifi
        fe
        control_time_nodes
        I_T
        I_u1
        I_v1
        I_v2
    end

    methods (Access = public)

        function [f] = f_AD(this, y, z, t)
            T = this.I_T * y;
            u1 = this.I_u1 * y;
            v1 = this.I_v1 * y;
            v2 = this.I_v2 * y;

            w = this.con_hifi.Temporal_Weights(t);
            zt = reshape(z, this.fe.m, this.con_hifi.control_time_nodes) * w;
            heating = this.fe.M * zt;

            T_coll = this.fe.nodes_to_coll_points * T;
            v1_coll = this.fe.nodes_to_coll_points * v1;
            v2_coll = this.fe.nodes_to_coll_points * v2;

            R_coll = this.con_hifi.Evaluate_Arrhenius_Law(T_coll);

            reaction_1 = this.con_hifi.react_rate_1 * this.fe.coll_point_integration * (R_coll .* v1_coll .* v2_coll);

            rhs_T = -this.con_hifi.diff_T * this.fe.S * T + heating - this.con_hifi.cooling_1 * reaction_1;
            rhs_u1 = -this.con_hifi.diff_u1 * this.fe.S * u1 + reaction_1;
            rhs_v1 = -this.con_hifi.diff_v1 * this.fe.S * v1 - reaction_1;
            rhs_v2 = -this.con_hifi.diff_v2 * this.fe.S * v2 - reaction_1;

            f_T = this.fe.M \ rhs_T;
            f_u1 = this.fe.M \ rhs_u1;
            f_v1 = this.fe.M \ rhs_v1;
            f_v2 = this.fe.M \ rhs_v2;

            f = this.I_T' * f_T + this.I_u1' * f_u1 + this.I_v1' * f_v1 + this.I_v2' * f_v2;
        end

        function [h] = h_AD(this, z)
            T_0 = this.con_hifi.Initial_T();
            u1_0 = this.con_hifi.Initial_u1();
            v1_0 = this.con_hifi.Initial_v1();
            v2_0 = this.con_hifi.Initial_v2();

            h = this.I_T' * T_0 + this.I_u1' * u1_0 + this.I_v1' * v1_0 + this.I_v2' * v2_0;
        end

    end

    methods (Access = public)

        function this = Thermochemical_LoFi_Constraint_AD(con_hifi)
            this@Dynamic_Constraint_AD(4 * con_hifi.fe.m, con_hifi.fe.m * con_hifi.control_time_nodes, con_hifi.T, con_hifi.n_t);

            this.con_hifi = con_hifi;
            this.fe = con_hifi.fe;
            this.control_time_nodes = con_hifi.control_time_nodes;

            this.I_T = zeros(this.fe.m, 4 * this.fe.m);
            this.I_T(:, 1:this.fe.m) = eye(this.fe.m);

            this.I_u1 = zeros(this.fe.m, 4 * this.fe.m);
            this.I_u1(:, (this.fe.m + 1):(2 * this.fe.m)) = eye(this.fe.m);

            this.I_v1 = zeros(this.fe.m, 4 * this.fe.m);
            this.I_v1(:, (2 * this.fe.m + 1):(3 * this.fe.m)) = eye(this.fe.m);

            this.I_v2 = zeros(this.fe.m, 4 * this.fe.m);
            this.I_v2(:, (3 * this.fe.m + 1):(4 * this.fe.m)) = eye(this.fe.m);
        end

    end

end
