classdef Adv_Diff < Constrained_ODE_Optimization


    properties
        A;
        M;
        x;
        beta_reg;
    end

    methods (Access = public)

         %% Instantiation of base class pure virtual functions for gradient computation
        function [val, grad_y] = Time_Instance_Objective(obj,y,t)
            target = obj.Evaluate_Target(t,obj.x);
            val = .5*(y-target)'*obj.M*(y-target);
            grad_y = obj.M*(y-target);
        end

        function [val,grad_z] = Regularization_Objective(obj,z)
            val = 0.0;
            grad_z = 0*z;
        end

        function [f, f_y, f_z] = Time_Instance_RHS(obj,y,z,t)
            % Extract the control for the given time.
            w = obj.Temporal_Weights(t);
            zt = reshape(z,obj.m,obj.N)*w;

            % Evaluate the RHS and its derivatives.
            f = linsolve(obj.M,-obj.A*y + obj.M*zt);
            f_y = linsolve(obj.M,-obj.A);
            f_z = kron(w',eye(obj.m));
        end

        function [h, h_z] = Initial_Condition(obj,z)
            h = zeros(obj.m,1);
            h_z = zeros(obj.m,obj.n);
        end

        %% Instantiation of base class pure virtual functions for hessian-vector product computation
        function [Mv] = Time_Instance_Objective_yy_Apply(obj,v,y,t)
            Mv = obj.M*v;
        end

        function [Mv] = Regularization_Objective_zz_Apply(obj,v,z)
            Mv = 0*v;
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
        function obj = Adv_Diff(m,n,T,Neumann)
            obj = obj@Constrained_ODE_Optimization(m,n,T,Neumann);

            Pe = 1; % Peclet number

            % Spatial domain
            obj.x = linspace(0,1,m)';
            h = obj.x(2)-obj.x(1);

            M = diag(4*ones(1,m)) + diag(ones(1,m-1),1) + diag(ones(1,m-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*h*M;

            S = diag(2*ones(1,m)) + (-1)*diag(ones(1,m-1),1) + (-1)*diag(ones(1,m-1),-1);
            S(1,1) = .5*S(1,1);
            S(end,end) = .5*S(end,end);
            S = (1/h)*S;

            V = diag(0*ones(1,m)) + (1/2)*diag(ones(1,m-1),1) + (-1/2)*diag(ones(1,m-1),-1);
            V(1,1) = -1/2;
            V(end,end) = 1/2;

            % Need to define A as the spatial discretization
            A = S + Pe*V;
            % f(y,z) = inv(M)*(-A*y + M*z)

            obj.A = A;
            obj.M = M;
            obj.beta_reg = 10^-5;
        end

        function [target] = Evaluate_Target(obj,t,x)
           target = t^2*exp(-50*(x-.5).^2);
        end

        function [w] = Temporal_Weights(obj,t)
            w = (obj.t_mesh-t)/(obj.t_mesh(2)-obj.t_mesh(1));
            Im = intersect(find(w<=0),find(abs(w)<=1));
            Ip = intersect(find(w>0),find(abs(w)<=1));
            I = find(abs(w)>1);
            w(I) = 0;
            w(Im) = 1+w(Im);
            w(Ip) = 1-w(Ip);
        end

    end

end
