classdef Adv_Diff_Gaussian_Source < Constrained_ODE_Optimization


    properties
        A;
        M;
        B;
        Br;
        x;
        num_space_control_nodes;
        time_weights;
        source_loc;
        z_time_mesh;
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
            val = .5*obj.beta_reg*z'*(kron(diag(obj.time_weights(2:end)),obj.Br'*obj.M*obj.Br)*z);
            grad_z = obj.beta_reg*(kron(diag(obj.time_weights(2:end)),obj.Br'*obj.M*obj.Br))*z;
        end

        function [f, f_y, f_z] = Time_Instance_RHS(obj,y,z,t)
            % Extract the control for the given time.
            w = obj.Temporal_Weights(t);
            zt = reshape(z,obj.num_space_control_nodes,length(obj.z_time_mesh))*w;

            % Evaluate the RHS and its derivatives.
            f = linsolve(obj.M,-obj.A*y + obj.M*obj.B*zt);
            f_y = linsolve(obj.M,-obj.A);
            f_z = kron(w',obj.B);
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
            Mv = obj.beta_reg*(kron(diag(obj.time_weights(2:end)),obj.Br'*obj.M*obj.Br))*v;
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
        function obj = Adv_Diff_Gaussian_Source(m,n,T,N,num_space_control_nodes)
            obj = obj@Constrained_ODE_Optimization(m,n,T,N);

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
            obj.beta_reg = 10^-3;

            % Trapazoid rule for time integration
            time_weights = ones(N,1);
            time_weights(1) = .5*time_weights(1);
            time_weights(end) = .5*time_weights(end);
            time_weights = T*time_weights/sum(time_weights);
            obj.time_weights = time_weights;

            % Control matrix (Gaussian)
            source_loc = linspace(0,1,num_space_control_nodes)';
            B = zeros(m,num_space_control_nodes);
            for k = 1:num_space_control_nodes
               B(:,k) = exp(-200*(obj.x-source_loc(k)).^2);
            end
            Br = B;
            B(1,:) = 0*B(1,:);
            B(end,:) = 0*B(end,:);

            obj.B = B;
            obj.Br = Br;
            obj.num_space_control_nodes = num_space_control_nodes;
            obj.source_loc = source_loc;
            obj.z_time_mesh = linspace(0,T,N)';
            obj.z_time_mesh = obj.z_time_mesh(2:end);
        end

        function [target] = Evaluate_Target(obj,t,x)
           target = 0.2*t^2*exp(-10*(x-.5).^2);
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
