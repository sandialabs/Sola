classdef Sparse_Matrix_Sqrt < Matrix_Sqrt

    properties
        A
    end

    %% Implementation of base class functions
    methods (Access = public)

        function [vec_out] = Matrix_Apply(this, vec_in)
            vec_out = this.A \ vec_in;
        end

    end

    %% Constructor
    methods

        function this = Sparse_Matrix_Sqrt(A)
            this.A = A;
        end

    end

end
