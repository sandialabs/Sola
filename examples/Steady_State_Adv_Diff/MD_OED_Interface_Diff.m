classdef MD_OED_Interface_Diff < MD_OED_Interface

    properties
        E
        M
        R
    end

    methods (Access = public)

        % Compute L*v, where Sigma = L*L^T
        function [L_v] = Apply_Design_Cov_Factor(this, v)
            L_v = linsolve(this.E, this.R' * v);
        end

        % Compute Sigma^{-1}*v
        function [Sigmainv_v] = Apply_Design_Cov_Inverse(this, v)
            tmp = this.E * v;
            tmp = linsolve(this.M, tmp);
            Sigmainv_v = this.E * tmp;
        end

    end

    methods

        function this = MD_OED_Interface_Diff(data_interface, con)
            this@MD_OED_Interface(data_interface);
            this.M = con.M;
            this.R = chol(this.M); % M = R^T*R
            this.E = (4) * ((1.e-2) * con.S + con.M); % Sigma = E^{-1}*M*E^{-1}
        end

    end

end
