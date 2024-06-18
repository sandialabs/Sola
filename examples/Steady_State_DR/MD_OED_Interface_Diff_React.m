classdef MD_OED_Interface_Diff_React < MD_OED_Interface

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

        function this = MD_OED_Interface_Diff_React(varargin)
            if nargin == 4
                [data_interface, con, alpha_zd, beta_zd] = deal(varargin{:});
            elseif nargin == 3
                [con, alpha_zd, beta_zd] = deal(varargin{:});
                data_interface = NaN;
            else
                error("Not enough (or too many) inputs.");
            end
            this@MD_OED_Interface(data_interface);
            this.M = con.M;
            this.R = chol(this.M); % M = R^T*R
            this.E = sqrt(alpha_zd) * (beta_zd * con.S + con.M); % Sigma = E^{-1}*M*E^{-1}
        end

    end

end
