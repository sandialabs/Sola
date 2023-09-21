classdef Diff_React_Opt < Constrained_Optimization
    
    
    properties
        diff_react_lofi;
        reg_coeff;
        m;
        M;
        T;
    end
    
    methods (Access = public)
        
        function [val, grad_u, grad_z] = Objective(obj,u,z)
            val = (1/2)*(u-obj.T)'*obj.M*(u-obj.T) + (1/2)*(obj.reg_coeff)*z'*obj.M*z;
            grad_u = obj.M*(u-obj.T);
            grad_z = (obj.reg_coeff)*obj.M*z;
        end
        
        function [u] = State_Solve(obj,z)
            u = obj.diff_react_lofi.State_Solve(z);
        end
        
        function [Mv] = c_u_Transpose_Inverse_Apply(obj,v,u,z)
            A = obj.diff_react_lofi.State_Jacobian();
            Mv = (A')\v;
        end
        
        function [Mv] = c_z_Transpose_Apply(obj,v,u,z)
            A = obj.diff_react_lofi.Control_Jacobian();
            Mv = A'*v;
        end
        
        function [Mv] = c_u_Inverse_Apply(obj,v,u,z)
            A = obj.diff_react_lofi.State_Jacobian();
            Mv = A\v;
        end
        
        function [Mv] = c_z_Apply(obj,v,u,z)
            A = obj.diff_react_lofi.Control_Jacobian();
            Mv = A*v;
        end
        
        function [Mv] = c_uu_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = c_uz_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = c_zu_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = c_zz_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = Objective_uu_Apply(obj,v,u,z)
            Mv = obj.M*v;
        end
        
        function [Mv] = Objective_uz_Apply(obj,v,u,z)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = Objective_zu_Apply(obj,v,u,z)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = Objective_zz_Apply(obj,v,u,z)
            Mv = obj.reg_coeff*obj.M*v;
        end
        
    end
    
    methods (Access = public)
         
        function obj = Diff_React_Opt(diff_react_lofi,reg_coeff)
            obj = obj@Constrained_Optimization();
            obj.diff_react_lofi = diff_react_lofi;
            obj.reg_coeff = reg_coeff;
            obj.m = size(obj.diff_react_lofi.A,1);
            obj.M = obj.diff_react_lofi.pde_meshing.M;
            x = obj.diff_react_lofi.pde_meshing.x;
            y = obj.diff_react_lofi.pde_meshing.y;
            obj.T = 15*exp(-5*( x.^2 + y.^2 )).*(1+.5*sin(pi*x).^2);
        end
        
    end
end

