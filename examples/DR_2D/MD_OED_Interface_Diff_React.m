classdef MD_OED_Interface_Diff_React < MD_OED_Interface

    properties
        E
        M
        R
    end

    methods (Access = public)

        % Compute L*v, where Sigma = L*L^T
        function [L_v] = Apply_Design_Cov_Factor(this, v)
            L_v = this.E \ (this.R' * v);
        end

        % Compute Sigma^{-1}*v
        function [Sigmainv_v] = Apply_Design_Cov_Inverse(this, v)
            tmp = this.E * v;
            tmp = this.M \ tmp;
            Sigmainv_v = this.E * tmp;
        end

    end

    methods

        function this = MD_OED_Interface_Diff_React(data_interface, con)
            this@MD_OED_Interface(data_interface);
            this.M = con.diff_react_lofi.M;
            this.R = chol(this.M); % M = R^T*R
            this.E = (3.e-2) * ((5.e-2) * con.diff_react_lofi.S + con.diff_react_lofi.M); % Sigma = E^{-1}*M*E^{-1}
        end

    end

end
