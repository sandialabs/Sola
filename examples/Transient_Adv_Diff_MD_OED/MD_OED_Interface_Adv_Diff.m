classdef MD_OED_Interface_Adv_Diff < MD_OED_Interface

    properties
        Sigmainv
        L
        scaling
    end

    methods (Access = public)

        % Compute L*v, where Sigma = L*L^T
        function [L_v] = Apply_Design_Cov_Factor(this, v)
            L_v = sqrt(this.scaling) * this.L * v;
        end

        % Compute Sigma^{-1}*v
        function [Sigmainv_v] = Apply_Design_Cov_Inverse(this, v)
            Sigmainv_v = (1 / this.scaling) * this.Sigmainv * v;
        end

    end

    methods

        function this = MD_OED_Interface_Adv_Diff(data_interface, obj)
            this@MD_OED_Interface(data_interface);

            this.Sigmainv = kron(diag(obj.time_weights(2:end)), obj.Br' * obj.M * obj.Br);
            R = chol(this.Sigmainv);
            this.L = linsolve(R, eye(size(R, 1)));

            this.scaling = 1.e-3;
        end

    end

end
