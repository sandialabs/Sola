classdef Diff_React_Constraint < Constraint
    
    
    properties
        diff_react_lofi;
        m;
    end
    
    methods (Access = public)
        
        function [u] = State_Solve(this,z)
            u = this.diff_react_lofi.State_Solve(z);
        end
        
        function [Mv] = c_u_Transpose_Inverse_Apply(this,v,u,z)
            A = this.diff_react_lofi.State_Jacobian();
            Mv = (A')\v;
        end
        
        function [Mv] = c_z_Transpose_Apply(this,v,u,z)
            A = this.diff_react_lofi.Control_Jacobian();
            Mv = A'*v;
        end
        
        function [Mv] = c_u_Inverse_Apply(this,v,u,z)
            A = this.diff_react_lofi.State_Jacobian();
            Mv = A\v;
        end
        
        function [Mv] = c_z_Apply(this,v,u,z)
            A = this.diff_react_lofi.Control_Jacobian();
            Mv = A*v;
        end
        
        function [Mv] = c_uu_Apply(this,v,u,z,lambda)
            Mv = zeros(this.m,1);
        end
        
        function [Mv] = c_uz_Apply(this,v,u,z,lambda)
            Mv = zeros(this.m,1);
        end
        
        function [Mv] = c_zu_Apply(this,v,u,z,lambda)
            Mv = zeros(this.m,1);
        end
        
        function [Mv] = c_zz_Apply(this,v,u,z,lambda)
            Mv = zeros(this.m,1);
        end

    end
    
    methods (Access = public)
         
        function this = Diff_React_Constraint(diff_react_lofi)
            this = this@Constraint();
            this.diff_react_lofi = diff_react_lofi;
            this.m = size(this.diff_react_lofi.A,1);
        end
        
    end
end

