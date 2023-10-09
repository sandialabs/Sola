classdef Example_5 < Constrained_ODE_Optimization

    % Solve the optimiation problem
    % min_{z} J(S(z),z) = int_0^T g(S(z)(t))dt + R(z)
    % where S(z) solves the ordinary differential equation
    % dy/dt = [2*t*y_1 ; 3*t^2*y_2 + t^2 - z(t)]
    % y(0) = [1, 1]
    % g(y) = (y_1-exp(t^2))^2 + (y_2-exp(t^3))^2
    % R(z) = int_0^T (z(t) - t^2)^2dt

    properties
        z_time_mesh;
        weights;
        beta_reg;
    end

    methods (Access = public)

         %% Instantiation of base class pure virtual functions for gradient computation
        function [val, grad_y] = Time_Instance_Objective(obj,y,t)
            val = (y(1)-exp(t^2))^2 + (y(2)-exp(t^3))^2;
            grad_y = zeros(2,1);
            grad_y(1) = 2*(y(1)-exp(t^2));
            grad_y(2) = 2*(y(2)-exp(t^3));
        end

        function [val,grad_z] = Regularization_Objective(obj,z)
            val = obj.beta_reg*(z-obj.z_time_mesh.^2)'*diag(obj.weights)*(z-obj.z_time_mesh.^2);
            grad_z = obj.beta_reg*2*diag(obj.weights)*(z-obj.z_time_mesh.^2);
        end

        function [f, f_y, f_z] = Time_Instance_RHS(obj,y,z,t)
            w = obj.Temporal_Weights(t);
            zt = w'*z;
            f = [ 2*t*y(1) ; 3*t^2*y(2) + t^2 - zt ];
            f_y = [ 2*t , 0 ; 0 , 3*t^2 ];
            f_z = zeros(2,obj.n);
            f_z(2,:) = -w';
        end

        function [h, h_z] = Initial_Condition(obj,z)
            h = [1 ; 1];
            h_z = zeros(2,length(z));
        end

        %% Instantiation of base class pure virtual functions for hessian-vector product computation
        function [Mv] = Time_Instance_Objective_yy_Apply(obj,v,y,t)
            A = 2*eye(2);
            Mv = A*v;
        end

        function [Mv] = Regularization_Objective_zz_Apply(obj,v,z)
            Mv = obj.beta_reg*2*diag(obj.weights)*v;
        end

        function [Mv] = Time_Instance_RHS_yy_Apply(obj,v,y,z,t,lambda)
            num_vecs = size(v,2);
            Mv = zeros(obj.m,num_vecs);
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(obj,v,y,z,t,lambda)
            num_vecs = size(v,2);
            Mv = zeros(obj.m,num_vecs);
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(obj,v,y,z,t,lambda)
            num_vecs = size(v,2);
            Mv = zeros(obj.n,num_vecs);
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(obj,v,y,z,t,lambda)
            num_vecs = size(v,2);
            Mv = zeros(obj.n,num_vecs);
        end

        function [Mv] = Initial_Condition_zz_Apply(obj,v,z,lambda)
            num_vecs = size(v,2);
            Mv = zeros(obj.n,num_vecs);
        end

    end

    methods (Access = public)
        function obj = Example_5(m,n,T,N)
            obj = obj@Constrained_ODE_Optimization(m,n,T,N);
            obj.z_time_mesh = linspace(0,T,n+1)';
            obj.z_time_mesh = obj.z_time_mesh(2:end);

            weights = ones(obj.n+1,1);
            weights(1) = .5; weights(end) = .5;
            weights = obj.T*weights/sum(weights);
            weights = weights(2:end);
            obj.weights = weights;
            obj.beta_reg = 10^-4;
        end

        function [w] = Temporal_Weights(obj,t)
           w = (obj.z_time_mesh-t)/(obj.z_time_mesh(2)-obj.z_time_mesh(1));
           Im = intersect(find(w<=0),find(abs(w)<=1));
           Ip = intersect(find(w>0),find(abs(w)<=1));
           I = find(abs(w)>1);
           w(I) = 0;
           w(Im) = 1+w(Im);
           w(Ip) = 1-w(Ip);
        end

    end

end
