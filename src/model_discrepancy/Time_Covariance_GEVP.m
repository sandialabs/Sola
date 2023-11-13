classdef Time_Covariance_GEVP < Randomized_GEVP

    properties
        E_tu
        E_td
        is_computed
    end

    methods (Access = public)

        function [vec_out] = Apply_Operator(this, vec_in)
            vec_out = this.E_tu \ vec_in;
        end

        function [vec_out] = Apply_Weighting_Operator(this, vec_in)
            vec_out = this.E_td \ vec_in;
        end

        function [vec_out] = Apply_Weighting_Operator_Inverse(this, vec_in)
            vec_out = this.E_td * vec_in;
        end

    end

    methods

        function this = Time_Covariance_GEVP(E_tu, E_td)
            vec = zeros(size(E_tu, 1), 1);
            this@Randomized_GEVP(vec);
            this.E_tu = E_tu;
            this.E_td = E_td;
        end

    end

end
