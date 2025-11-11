classdef Sparse_Matrix_Sqrt < Matrix_Sqrt

    properties
        A
        L
    end

    %% Implementation of base class functions
    methods (Access = public)

        function [vec_out] = Matrix_Apply(this, vec_in)
            vec_out = this.A * vec_in;
        end

        function [vec_out] = Preconditioner_Apply(this, vec_in)
            vec_out = this.L \ vec_in;
        end

        function [vec_out] = Preconditioner_Transpose_Apply(this, vec_in)
            vec_out = this.L' \ vec_in;
        end

        function [vec_out] = Preconditioner_Inverse_Apply(this, vec_in)
            vec_out = this.L * vec_in;
        end

    end

    %% Constructor
    methods

        function this = Sparse_Matrix_Sqrt(A,L)
            arguments
                A {mustBeNumeric}
                L {mustBeNumeric} = 1;
            end
            this.A = A;
            this.L = L;
        end

    end

end
