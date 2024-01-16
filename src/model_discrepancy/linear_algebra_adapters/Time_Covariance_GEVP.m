classdef Time_Covariance_GEVP < Randomized_GEVP

    properties
        E_t
        M_t
        is_computed
    end

    methods (Access = public)

        function [vec_out] = Apply_Operator(this, vec_in)
            vec_out = this.E_t \ vec_in;
        end

        function [vec_out] = Apply_Weighting_Operator(this, vec_in)
            vec_out = this.M_t \ vec_in;
        end

        function [vec_out] = Apply_Weighting_Operator_Inverse(this, vec_in)
            vec_out = this.M_t * vec_in;
        end

    end

    methods

        function this = Time_Covariance_GEVP(E_t, M_t)
            vec = zeros(size(E_t, 1), 1);
            this@Randomized_GEVP(vec);
            this.E_t = E_t;
            this.M_t = M_t;
        end

    end

end
