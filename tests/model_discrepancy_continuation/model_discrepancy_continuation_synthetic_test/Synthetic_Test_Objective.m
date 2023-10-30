classdef Synthetic_Test_Objective < Objective
    
    properties
        con;
        M;
        R;
        uT;
        m;
        n;
    end
    
    methods (Access = public)
        
        function [val, grad_u, grad_z] = J(this, u, z)
            val = (1/2)*(u-this.uT)'*this.M*(u-this.uT) + (1/2)*z'*this.R*z;
            grad_u = this.M*(u-this.uT);
            grad_z = this.R*z;
        end
        
        function [Mv] = J_uu_Apply(this, v, u, z)
            Mv = this.M*v;
        end
        
        function [Mv] = J_uz_Apply(this, v, u, z)
            Mv = zeros(this.m,size(v,2));
        end
        
        function [Mv] = J_zu_Apply(this, v, u, z)
            Mv = zeros(this.n,size(v,2));
        end
        
        function [Mv] = J_zz_Apply(this, v, u, z)
            Mv = this.R*v;
        end
        
    end
    
    methods (Access = public)
        
        function this = Synthetic_Test_Objective(con)
            this.con = con;
            m = con.m;
            n = con.n;
            M = 2*eye(m);
            R = (1.e-2)*eye(n);
            zT = linspace(0,1,n)';
            uT = con.State_Solve(zT);
            this.m = m;
            this.n = n;
            this.M = M;
            this.R = R;
            this.uT = uT;
        end
        
    end
end
