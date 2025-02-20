classdef Matrix_Sqrt_Test < Matrix_Sqrt

    properties
        M
    end

    methods (Access = public)

        function [vec_out] = Matrix_Apply(this, vec_in)
            vec_out = this.M * vec_in;
        end

    end

    methods

        function this = Matrix_Sqrt_Test(m)
            h = 1 / (m - 1);
            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;
        end

    end

end
