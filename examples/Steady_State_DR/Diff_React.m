classdef Diff_React < Constrained_Optimization
    
    
    properties
        m;
        reg_coeff;
        reg_mat;
        T;
        diff_coeff;
        react_coeff;
        x;
        M;
        S;
        coll_weights;
        coll_points;
        nodes_to_coll_points;
    end
    
    methods (Access = public)
        
        %% Pure virtual functions for gradient computation
        
        function [val, grad_u, grad_z] = Objective(obj,u,z)
            val = (1/2)*(u-obj.T)'*obj.M*(u-obj.T) + (1/2)*(obj.reg_coeff)*z'*obj.reg_mat*z;
            grad_u = obj.M*(u-obj.T);
            grad_z = (obj.reg_coeff)*obj.reg_mat*z;
        end
        
        function [u] = State_Solve(obj,z)
            % Solve system without reaction term to generate nonlinear solver initial iterate
            A = obj.diff_coeff*obj.S + (-obj.react_coeff)*(1.e-3)*eye(size(obj.S,1));
            b = obj.M*z;
            u0 = linsolve(A,b);
            
            % Execute nonlinear solve to determine the state
            options = optimoptions('fsolve','Display','none','OptimalityTolerance',1.e-14,'SpecifyObjectiveGradient',true,'CheckGradients',false);
            u = fsolve(@(u)obj.Constraint_Evaluation(u,z),u0,options);
        end
        
        function [c,c_u,c_z] = Constraint_Evaluation(obj,u,z)
            [R,R_u] = obj.Assemble_Reaction_Function(u);
            c =  obj.diff_coeff*obj.S*u - obj.react_coeff*R - obj.M*z;
            c_u = obj.diff_coeff*obj.S - obj.react_coeff*R_u;
            c_z = -obj.M;
        end
        
        function [c_uu] = Constraint_Hessian(obj,u,lambda)
            c_uu = -obj.react_coeff*obj.Assemble_Reaction_Function_Hessian(u,lambda);
        end
        
        function [R,R_u] = Assemble_Reaction_Function(obj,u)
            u_nodes = obj.nodes_to_coll_points*u;
            [R_nodes,R_prime_nodes] = obj.Reaction_Function(u_nodes);
            R = obj.nodes_to_coll_points'*(obj.coll_weights.*R_nodes);
            R_u = obj.nodes_to_coll_points'*(diag(obj.coll_weights)*R_prime_nodes)*obj.nodes_to_coll_points;
        end
        
        function [R_uu] = Assemble_Reaction_Function_Hessian(obj,u,lambda)
            u_nodes = obj.nodes_to_coll_points*u;
            lambda_nodes = obj.nodes_to_coll_points*lambda;
            R_uu = obj.nodes_to_coll_points'*diag(obj.coll_weights)*obj.Reaction_Function_Hessian(u_nodes,lambda_nodes)*obj.nodes_to_coll_points;
        end
        
        function [R,R_prime] = Reaction_Function(obj,u,x)
            R = u.^2;
            R_prime = 2*diag(u);
        end
        
        function [R_prime_prime] = Reaction_Function_Hessian(obj,u,lambda,x)
            R_prime_prime = 2*diag(lambda);
        end
        
        function [diff] = Finite_Difference_Reaction_Function_Jacobian(obj,u)
            [R,R_u] = obj.Assemble_Reaction_Function(u);
            h = 10.^(-1:-1:-6);
            v = randn(length(u),1);
            v = v/norm(v);
            diff = zeros(6,1);
            for k = 1:6
                [R_pert,~] = obj.Assemble_Reaction_Function(u+h(k)*v);
                diff(k) = norm(R_u*v - (R_pert - R)/h(k))/norm(R_u*v);
            end
            disp(log10(diff'))
        end
        
        function [diff] = Finite_Difference_Reaction_Function_Hessian(obj,u,lambda)
            R_uu = obj.Assemble_Reaction_Function_Hessian(u,lambda);
            [~,R_u] = obj.Assemble_Reaction_Function(u);
            h = 10.^(-1:-1:-6);
            v = randn(length(u),1);
            v = v/norm(v);
            diff = zeros(6,1);
            for k = 1:6
                [~,R_u_pert] = obj.Assemble_Reaction_Function(u+h(k)*v);
                diff(k) = norm(R_uu*v - (R_u_pert'*lambda - R_u'*lambda)/h(k))/norm(R_uu*v);
            end
            disp(log10(diff'))
        end
        
        function [diff] = Finite_Difference_Constraint_Hessian(obj,u,z,lambda)
            c_uu = obj.Constraint_Hessian(u,lambda);
            [~,c_u,~] = obj.Constraint_Evaluation(u,z);
            h = 10.^(-1:-1:-6);
            v = randn(length(u),1);
            v = v/norm(v);
            diff = zeros(6,1);
            for k = 1:6
                [~,c_u_pert,~] = obj.Constraint_Evaluation(u+h(k)*v,z);
                diff(k) = norm(c_uu*v - (c_u_pert'*lambda - c_u'*lambda)/h(k))/norm(c_uu*v);
            end
            disp(log10(diff'))
        end
        
        function [Mv] = c_u_Transpose_Inverse_Apply(obj,v,u,z)
            [~,c_u] = obj.Constraint_Evaluation(u,z);
            Mv = linsolve(c_u',v);
        end
        
        function [Mv] = c_z_Transpose_Apply(obj,v,u,z)
            [~,~,c_z] = obj.Constraint_Evaluation(u,z);
            Mv = c_z'*v;
        end
        
        function [Mv] = c_u_Inverse_Apply(obj,v,u,z)
            [~,c_u] = obj.Constraint_Evaluation(u,z);
            Mv = linsolve(c_u,v);
        end
        
        function [Mv] = c_z_Apply(obj,v,u,z)
            [~,~,c_z] = obj.Constraint_Evaluation(u,z);
            Mv = c_z*v;
        end
        
        function [Mv] = c_uu_Apply(obj,v,u,z,lambda)
            c_uu = obj.Constraint_Hessian(u,lambda);
            Mv = c_uu*v;
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
         
        function obj = Diff_React(m,diff_coeff,react_coeff,reg_coeff)
            obj = obj@Constrained_Optimization();
            obj.m = m;
            obj.diff_coeff = diff_coeff;
            obj.react_coeff = react_coeff;
            obj.x = linspace(0,1,m)';
            obj.reg_coeff = reg_coeff;
            obj.T = 20*(obj.x+.5).*(1.3-obj.x); %6*(obj.x.*(1-obj.x));
            
            h = obj.x(2)-obj.x(1);
            
            coll_points = zeros(2*(m-1),1);
            for k = 1:(m-1)
                map_to_coll = (1:2)' + 2*(k-1);
               coll_points(map_to_coll) = obj.x(k) + h*((1/sqrt(3))*[-1;1] + 1)/2;  
            end
            coll_weights = (h/2)*ones(2*(m-1),1);
            
            nodes_to_coll_points = zeros(2*(m-1),m);
            for k = 1:(m-1)
                map_to_coll = (1:2)' + 2*(k-1);
               nodes_to_coll_points(map_to_coll(1),k) = (coll_points(map_to_coll(1))-obj.x(k+1))/(obj.x(k)-obj.x(k+1));
               nodes_to_coll_points(map_to_coll(1),k+1) = (coll_points(map_to_coll(1))-obj.x(k))/(obj.x(k+1)-obj.x(k));
               nodes_to_coll_points(map_to_coll(2),k) = (coll_points(map_to_coll(2))-obj.x(k+1))/(obj.x(k)-obj.x(k+1));
               nodes_to_coll_points(map_to_coll(2),k+1) = (coll_points(map_to_coll(2))-obj.x(k))/(obj.x(k+1)-obj.x(k));
            end
            
            obj.coll_weights = coll_weights;
            obj.coll_points = coll_points;
            obj.nodes_to_coll_points = nodes_to_coll_points;
            
            M = diag(4*ones(1,m)) + diag(ones(1,m-1),1) + diag(ones(1,m-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*h*M;
            obj.M = M;
            
            S = diag(2*ones(1,m)) + (-1)*diag(ones(1,m-1),1) + (-1)*diag(ones(1,m-1),-1);
            S(1,1) = .5*S(1,1);
            S(end,end) = .5*S(end,end);
            S = (1/h)*S;
            obj.S = S;
            
            obj.reg_mat = M;
        end
        
    end
end

