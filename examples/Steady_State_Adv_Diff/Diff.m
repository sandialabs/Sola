classdef Diff < Constrained_Optimization
    
    
    properties
        m;
        reg_coeff;
        reg_mat;
        T;
        diff_coeff;
        robin_coeff;
        x;
        M;
        S;
        robin_bc;
    end
    
    methods (Access = public)
        
        %% Pure virtual functions for gradient computation
        
        function [val, grad_u, grad_z] = Objective(obj,u,z)
            val = (1/2)*(u-obj.T)'*obj.M*(u-obj.T) + (1/2)*(obj.reg_coeff)*z'*obj.reg_mat*z;
            grad_u = obj.M*(u-obj.T);
            grad_z = (obj.reg_coeff)*obj.reg_mat*z;
        end
        
        function [u] = State_Solve(obj,z)
            A = obj.diff_coeff*obj.S + obj.robin_coeff*obj.robin_bc;
            b = (10^2)*obj.M*z;
            u = linsolve(A,b);
        end
        
        function [Mv] = c_u_Transpose_Inverse_Apply(obj,v,u,z)
            A = obj.diff_coeff*obj.S + obj.robin_coeff*obj.robin_bc;
            Mv = linsolve(A',v);
        end
        
        function [Mv] = c_z_Transpose_Apply(obj,v,u,z)
            Mv = -(10^2)*obj.M'*v;
        end
        
        function [Mv] = c_u_Inverse_Apply(obj,v,u,z)
            A = obj.diff_coeff*obj.S + obj.robin_coeff*obj.robin_bc;
            Mv = linsolve(A,v);
        end
        
        function [Mv] = c_z_Apply(obj,v,u,z)
            Mv = -(10^2)*obj.M*v;
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
            Mv = obj.reg_coeff*obj.reg_mat*v;
        end
        
    end
    
    methods (Access = public)
        
        function obj = Diff(adv_diff)
            obj = obj@Constrained_Optimization();
            obj.m = adv_diff.m;
            obj.x = adv_diff.x;
            obj.T = adv_diff.T;
            obj.reg_coeff = adv_diff.reg_coeff;
            obj.diff_coeff = adv_diff.diff_coeff;
            obj.robin_coeff = adv_diff.robin_coeff;
            obj.M = adv_diff.M;
            obj.S = adv_diff.S;
            obj.robin_bc = adv_diff.robin_bc;
            obj.reg_mat = adv_diff.reg_mat;
        end
        
    end
end

