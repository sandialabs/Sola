classdef Example_2 < Constrained_ODE_Optimization
    % Solve the optimiation problem
    % min_{z} J(S(z),z) = int_0^T g(S(z)(t))dt + R(z)
    % where S(z) solves the odinary differential equation
    % dy/dt = [y_1 , y_2]
    % y(0) = [z_1 , z_2]
    % g(y) = (y_1-exp(t))^2 + (y_2-exp(t))^2
    % R(z) = (z_1-1)^2 + (z_2-1)^2

    properties

    end

    methods (Access = public)

        %% Instantiation of base class pure virtual functions for gradient computation
        function [val, grad_y] = Time_Instance_Objective(obj,y,t)
            val = (y(1)-exp(t))^2 + (y(2)-exp(t))^2;
            grad_y = 2*(y-exp(t));
        end

        function [val,grad_z] = Regularization_Objective(obj,z)
            val = (z(1)-1)^2 + (z(2)-1)^2;
            grad_z = 2*(z-1);
        end

        function [f, f_y, f_z] = Time_Instance_RHS(obj,y,z,t)
            f = y;
            f_y = eye(2);
            f_z = zeros(2,2);
        end

        function [h, h_z] = Initial_Condition(obj,z)
            h = z;
            h_z = eye(2);
        end

        %% Instantiation of base class pure virtual functions for hessian-vector product computation
        function [Mv] = Time_Instance_Objective_yy_Apply(obj,v,y,t)
           A = 2*eye(2);
           Mv = A*v;
        end

        function [Mv] = Regularization_Objective_zz_Apply(obj,v,z)
            A = 2*eye(2);
            Mv = A*v;
        end

        function [Mv] = Time_Instance_RHS_yy_Apply(obj,v,y,z,t,lambda)
            Mv = 0*v;
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(obj,v,y,z,t,lambda)
            num_vecs = size(v,2);
            Mv = zeros(obj.m,num_vecs);
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(obj,v,y,z,t,lambda)
            num_vecs = size(v,2);
            Mv = zeros(length(z),num_vecs);
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(obj,v,y,z,t,lambda)
            Mv = 0*v;
        end

        function [Mv] = Initial_Condition_zz_Apply(obj,v,z,lambda)
            Mv = 0*z;
        end

    end

    methods (Access = public)
        function obj = Example_2(m,n,T,N)
            obj = obj@Constrained_ODE_Optimization(m,n,T,N);
        end

    end

end
