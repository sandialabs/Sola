classdef Synthetic_Test_Constraint < Constraint
    
    properties
        A;
        B;
        m;
        n;
    end
    
    methods (Access = public)
        
        function [u] = State_Solve(this, z)
            u = linsolve(this.A,this.B*z);
        end
        
        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            Mv = linsolve(this.A',v);
        end
        
        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            Mv = -this.B'*v;
        end
        
        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            Mv = linsolve(this.A,v);
        end
        
        function [Mv] = c_z_Apply(this, v, u, z)
            Mv = -this.B*v;
        end
        
        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m,size(v,2));
        end
        
        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m,size(v,2));
        end
        
        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.n,size(v,2));
        end
        
        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.n,size(v,2));
        end
        
    end
    
    methods (Access = public)

        function this = Synthetic_Test_Constraint( )
            m = 10;
            n = 9;
            
            A = eye(m);
            A(1,2) = -1;
            A(3,1) = -1;
            A(8,7) = -1;
            
            B = diag(1:m);
            B = B(:,1:n);
            B(m,n) = 10;
            
            this.m = m;
            this.n = n;
            this.A = A;
            this.B = B;
        end

    end
end
