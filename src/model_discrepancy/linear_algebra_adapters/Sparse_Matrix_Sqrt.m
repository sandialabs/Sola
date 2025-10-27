classdef Sparse_Matrix_Sqrt < Matrix_Sqrt

    properties
        A
        use_inverse
    end

    %% Implementation of base class functions
    methods (Access = public)

        function [vec_out] = Matrix_Apply(this, vec_in)
            if this.use_inverse
                vec_out = this.A \ vec_in;
            else
                vec_out = this.A * vec_in;
            end
        end

    end

    %% Constructor
    methods

        function this = Sparse_Matrix_Sqrt(A,use_inverse)
            this.A = A;
            this.use_inverse = use_inverse;
        end

    end

end
